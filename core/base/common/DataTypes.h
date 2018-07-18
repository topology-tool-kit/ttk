/// \ingroup base
/// \class ttk::DataTypes
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date May 2018.
///
///\brief TTK base package defining the standard types.

#ifndef _DATATYPES_H
#define _DATATYPES_H

namespace ttk{
  /// \brief Identifier type for simplices of any dimension.
  using LongSimplexId = long long int;

  /// \brief Identifier type for simplices of any dimension.
#ifdef TTK_USE_64BIT_IDS
  using SimplexId = long long int;
#else
  using SimplexId = int;
#endif

  /// \brief Identifier type for threads (i.e. with OpenMP).
  using ThreadId = int;

  /// \brief Identifier type for tasks (i.e. with OpenMP).
  using TaskId = int;

  /// default name for offset scalar field
  const char OffsetScalarFieldName[]="ttkOffsetScalarField";
}

#endif // _DATATYPES_H
