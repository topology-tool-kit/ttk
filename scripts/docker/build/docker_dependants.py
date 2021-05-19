#!/usr/bin/python3
#
# usage: python3 docker_descendants.py <image_id> ...

import sys
from subprocess import check_output


def main(images):
    image_ids = set(images)
    all_images = docker_images('--all', '--quiet')
    all_links = parse_links(docker_links(all_images))
    descendants = desc(image_ids, all_links)
    pred = lambda s: lambda line: s[:12] in line
    match = list(map(pred, descendants))
    return filter(lambda i: any(s(i) for s in match), docker_images())


def parse_links(lines):
    parseid = lambda s: s.replace('sha256:', '')
    for line in reversed(list(lines)):
        yield list(map(parseid, line.split()))


def desc(image_ids, links):
    if links:
        link, *tail = links
        if len(link) > 1:
            image_id, parent_id = link
            checkid = lambda i: parent_id.startswith(i)
            if any(map(checkid, image_ids)):
                return desc(image_ids | {image_id}, tail)
        return desc(image_ids, tail)
    return image_ids


def docker_links(images):
    cmd = [ 'docker', 'inspect', '--format={{.Id}} {{.Parent}}']
    return run(cmd + images)


def docker_images(*args):
    return run(('docker', 'images') + args)


def run(cmd):
    return check_output(cmd, universal_newlines=True).splitlines()


if __name__ == '__main__':
    print('\n'.join(main(sys.argv[1:])))


# License
# --------
# 
# Copyright (C) 2017 Aryeh Leib Taurog <python@aryehleib.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
