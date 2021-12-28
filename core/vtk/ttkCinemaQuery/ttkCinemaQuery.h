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
/// VTK wrapping code for the ttk::CinemaQuery package.
///
/// \param Input Input table (vtkTable)
/// \param Output Output table (vtkTable)
///
/// sa ttk::CinemaQuery
///
/// \b Online \b examples: \n
///   - <a href="https://topology-tool-kit.github.io/examples/cinemaIO/">Cinema
///   IO example</a> \n

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

  vtkSetMacro(SQLStatement, const std::string &);
  vtkGetMacro(SQLStatement, std::string);

  vtkSetMacro(ExcludeColumnsWithRegexp, bool);
  vtkGetMacro(ExcludeColumnsWithRegexp, bool);

  vtkSetMacro(RegexpString, const std::string &);
  vtkGetMacro(RegexpString, std::string);

protected:
  ttkCinemaQuery();
  ~ttkCinemaQuery() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::string SQLStatement{"SELECT * FROM InputTable0"};
  bool ExcludeColumnsWithRegexp{false};
  std::string RegexpString{".*"};
};
