#!/usr/bin/env python3
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

CLASS_PREFIX = NAME_PREFIX = 'ttk'
LABEL_PREFIX = 'TTK '

PARAVIEW_XML_DIR = '../paraview'
PARAVIEW_XML_OUTPUT = 'replace'  # 'replace', 'reformat', 'none'

COMPATIBILITY_XML_FILE = '../paraview/Compatibility/Compatibility.xml'
COMPATIBILITY_DEFINITION = '''
        <SourceProxy base_proxygroup="{group}"
                     base_proxyname="{proxy_name}"
                     class="{proxy_class}"
                     name="{proxy_name_old}"
                     label="{label} (deprecated)">
            <Deprecated deprecated_in="5.8" to_remove_in="5.9">
                Please update your state file to use {proxy_name} instead of {proxy_name_old}.
            </Deprecated>
        </SourceProxy> 
'''

STATE_FILE_DIR = '../../ttk-data/states'
STATE_FILE_OUTPUT = 'replace'  # 'replace', 'reformat', 'none'


def indent(elem, indentation=4 * ' ', level=0):
    i = '\n' + level * indentation
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + indentation
        for e in elem:
            indent(e, indentation, level + 1)
            if not e.tail or not e.tail.strip():
                e.tail = i + indentation
        if not e.tail or not e.tail.strip():
            e.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
    if len(elem.attrib) > 3:
        elem.attrib = {(i + indentation + k): v for k, v in elem.attrib.items()}


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.abspath(os.path.join(script_dir, PARAVIEW_XML_DIR))

    xml_files = Path(base_dir).rglob('*.xml')

    filters_map = {}
    compatiblity_root = ET.Element('ServerManagerConfiguration')

    if COMPATIBILITY_XML_FILE:
        compatiblity_file = os.path.realpath(os.path.join(script_dir, COMPATIBILITY_XML_FILE))
        try:
            e = ET.parse(compatiblity_file)
            compatiblity_root = e.getroot()
        except (ET.ParseError, FileNotFoundError):
            pass
    else:
        compatiblity_file = None

    for xml_file in xml_files:
        if os.path.realpath(xml_file) == compatiblity_file:
            continue
        try:
            e = ET.parse(xml_file)
        except ET.ParseError:
            continue
        if PARAVIEW_XML_OUTPUT == 'replace':
            with open(xml_file, 'r') as f:
                contents = f.read()

        for group in e.getroot().iter('ProxyGroup'):
            group_name = group.get('name')
            compatibility_group = compatiblity_root.find('./ProxyGroup[@name="{}"]'.format(group_name))
            if compatibility_group is None:
                compatibility_group = ET.SubElement(compatiblity_root, 'ProxyGroup', {'name': group_name})
            for proxy in group:
                proxy_name = proxy.get('name')
                proxy_class = proxy.get('class')
                proxy_label = proxy.get('label')

                if not proxy_class.startswith(CLASS_PREFIX):
                    print('Unexpected class name: ' + proxy_class)

                if proxy_label is not None and not proxy_label.startswith(LABEL_PREFIX):
                    proxy_label_old = proxy_label[:]
                    proxy_label = LABEL_PREFIX + proxy_label

                    if PARAVIEW_XML_OUTPUT == 'reformat':
                        proxy.set('label', proxy_label)
                    elif PARAVIEW_XML_OUTPUT == 'replace':
                        contents = contents.replace('label="{}"'.format(proxy_label_old), 'label="{}"'.format(proxy_label))

                    print('{}: Relabeled "{}" to "{}".'.format(xml_file.name, proxy_label_old, proxy_label))

                if not proxy_name.startswith(NAME_PREFIX):
                    proxy_name_old = proxy_name[:]
                    proxy_name = NAME_PREFIX + proxy_name

                    if PARAVIEW_XML_OUTPUT == 'reformat':
                        proxy.set('name', proxy_name)
                    elif PARAVIEW_XML_OUTPUT == 'replace':
                        contents = contents.replace('name="{}"'.format(proxy_name_old), 'name="{}"'.format(proxy_name))

                    compatibility_definition = COMPATIBILITY_DEFINITION.format(
                        group=group_name, proxy_class=proxy_class, proxy_name=proxy_name, proxy_name_old=proxy_name_old, label=proxy_label)

                    compatibility_group.append(ET.fromstring(compatibility_definition))

                    print('{}: Renamed {} to {}.'.format(xml_file.name, proxy_name_old, proxy_name))

                filters_map[proxy_name[len(NAME_PREFIX):]] = proxy_name
        if PARAVIEW_XML_OUTPUT == 'reformat':
            indent(e.getroot())
            e.write(xml_file)
        elif PARAVIEW_XML_OUTPUT == 'replace':
            with open(xml_file, 'w') as f:
                f.write(contents)

    if COMPATIBILITY_XML_FILE:
        indent(compatiblity_root)
        compatiblity = ET.ElementTree(compatiblity_root)
        compatiblity.write(os.path.join(script_dir, COMPATIBILITY_XML_FILE))


    if len(sys.argv) > 1:
        state_files = sys.argv[1:]
        state_files = [Path(f) for f in state_files]
    else:
        state_dir = os.path.abspath(os.path.join(script_dir, STATE_FILE_DIR))
        state_files = Path(state_dir).rglob('*.pvsm')

    for state_file in state_files:
        e = ET.parse(state_file)
        if STATE_FILE_OUTPUT == 'replace':
            with open(state_file, 'r') as f:
                contents = f.read()

        modified = False
        for proxy in e.getroot().iter('Proxy'):
            proxy_type = proxy.get('type')
            if proxy_type in filters_map:
                if STATE_FILE_OUTPUT == 'reformat':
                    proxy.set('type', filters_map[proxy_type])
                elif STATE_FILE_OUTPUT == 'replace':
                    contents = contents.replace('type="{}"'.format(proxy_type), 'type="{}"'.format(filters_map[proxy_type]))
                modified = True
                print('{}: Renamed "{}" to "{}".'.format(state_file.name, proxy_type, filters_map[proxy_type]))

        if modified:
            if STATE_FILE_OUTPUT in ['reformat', 'replace']:
                backup = state_file.with_suffix('.bak')
                if backup.exists():
                    print('Backup file {} already exists. Not writing file {}.'.format(backup, state_file))
                    continue
                state_file.rename(backup)

            if STATE_FILE_OUTPUT == 'reformat':
                e.write(state_file)
            elif STATE_FILE_OUTPUT == 'replace':
                with open(state_file, 'w') as f:
                    f.write(contents)


if __name__ == '__main__':
    main()
