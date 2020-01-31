/// \ingroup vtk
/// \class ttkCinemaQuery
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.09.2018
///
/// \brief TTK VTK-filter that uses a SQL statement to select a subset of a
/// vtkTable.
///
/// This filter creates a temporary SQLite3 database from the input table,
/// performs a SQL query, and then returns the result as a vtkTable.
///
/// VTK wrapping code for the @CinemaQuery package.
///
/// \param Input Input table (vtkTable)
/// \param Output Output table (vtkTable)
///
/// sa ttk::CinemaQuery

#pragma once

// VTK Module
#include <ttkCinemaQueryModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <CinemaQuery.h>

class TTKCINEMAQUERY_EXPORT ttkCinemaQuery : public ttkAlgorithm,
                                             protected ttk::CinemaQuery {
public:
  static ttkCinemaQuery *New();
  vtkTypeMacro(ttkCinemaQuery, ttkAlgorithm);

  vtkSetMacro(SQLStatement, std::string);
  vtkGetMacro(SQLStatement, std::string);

protected:
  ttkCinemaQuery();
  ~ttkCinemaQuery();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string SQLStatement{"SELECT * FROM InputTable0"};
};
