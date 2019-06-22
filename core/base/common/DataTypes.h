/// \ingroup base
/// \class ttk::DataTypes
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date May 2018.
///
///\brief TTK base package defining the standard types.

#ifndef _DATATYPES_H
#define _DATATYPES_H

namespace ttk {
  /// \brief Identifier type for simplices of any dimension.
  using LongSimplexId = long long int;

  /// \brief Identifier type for simplices of any dimension.
#ifdef TTK_ENABLE_64BIT_IDS
  using SimplexId = long long int;
#else
  using SimplexId = int;
#endif

  /// \brief Identifier type for threads (i.e. with OpenMP).
  using ThreadId = int;

  /// \brief Identifier type for tasks (i.e. with OpenMP).
  using TaskId = int;

  /// default name for mask scalar field
  const char MaskScalarFieldName[] = "ttkMaskScalarField";

  /// default name for vertex scalar field
  const char VertexScalarFieldName[] = "ttkVertexScalarField";

  /// default name for offset scalar field
  const char OffsetScalarFieldName[] = "ttkOffsetScalarField";

  /// default name for bivariate offset fields
  const char OffsetFieldUName[] = "ttkOffsetFieldU";
  const char OffsetFieldVName[] = "ttkOffsetFieldV";

  /// default value for critical index
  enum class CriticalType {
    Local_minimum = 0,
    Saddle1,
    Saddle2,
    Local_maximum,
    Degenerate,
    Regular
  };
  /// number of different critical types
  const int CriticalTypeNumber = 6;

} // namespace ttk

#endif // _DATATYPES_H
