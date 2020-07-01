/// \ingroup vtk
/// \class ttkContourTreeAlignment
/// \author Florian Wetzels (f_wetzels13@cs.uni-kl.de), Jonas Lukasczyk (jl@jluk.de)
/// \date 28.01.2020
///
/// \brief TTK VTK-filter that computes an alignment for a multiblock of contourtrees
///
/// VTK wrapping code for the @ContourTreeAlignment package.
///
/// TODO
///
/// \param Input TODO
/// \param Output TODO
///
/// \sa ttk::ContourTreeAlignment

#pragma once

// VTK Module
#include <ttkContourTreeAlignmentModule.h>

// VTK includes
#include <ttkAlgorithm.h>
#include <vtkInformation.h>

// TTK includes
#include <ContourTreeAlignment.h>

class TTKCONTOURTREEALIGNMENT_EXPORT ttkContourTreeAlignment
: public ttkAlgorithm,
protected ttk::ContourTreeAlignment{

public:

    // VTK stuff
    static ttkContourTreeAlignment* New();
    vtkTypeMacro(ttkContourTreeAlignment, ttkAlgorithm);

    // auto generated setters and getters
    vtkSetMacro(RandomSeed, int);
    vtkGetMacro(RandomSeed, int);
    vtkSetMacro(MatchTime, bool);
    vtkGetMacro(MatchTime, bool);
    vtkSetMacro(ArcMatchMode, int);
    vtkGetMacro(ArcMatchMode, int);
    vtkSetMacro(AlignmenttreeType, int);
    vtkGetMacro(AlignmenttreeType, int);
    vtkSetMacro(WeightCombinatorialMatch, float);
    vtkGetMacro(WeightCombinatorialMatch, float);
    vtkSetMacro(WeightArcMatch, float);
    vtkGetMacro(WeightArcMatch, float);
    vtkSetMacro(WeightScalarValueMatch, float);
    vtkGetMacro(WeightScalarValueMatch, float);
    vtkSetMacro(ExportJSON, bool);
    vtkGetMacro(ExportJSON, bool);
    vtkSetMacro(ExportPath, std::string);
    vtkGetMacro(ExportPath, std::string);

protected:
    // Specify the input data type of each input port
    int FillInputPortInformation(int port, vtkInformation *info) override;

    // Specify the data object type of each output port
    int FillOutputPortInformation(int port, vtkInformation *info) override;

    // Pass VTK data to the base code and convert base code output to VTK
    int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;

    // filter constructor and destructor
    ttkContourTreeAlignment();
    ~ttkContourTreeAlignment() override {};

private:
    // filter parameters
    int RandomSeed{1};
    bool MatchTime{false};
    int AlignmenttreeType{0};
    bool ExportJSON{false};
    std::string ExportPath{""};
    int ArcMatchMode{2};
    float WeightCombinatorialMatch{0};
    float WeightArcMatch{1};
    float WeightScalarValueMatch{0};

};
