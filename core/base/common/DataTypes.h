/// \ingroup base
/// \class ttk::DataTypes
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date May 2018.
///
///\brief TTK base package defining the standard types.

#ifndef DATATYPES_H
#define DATATYPES_H

namespace ttk{
  /// \brief Identifier type for simplices of any dimension.
  using SimplexId = long long int;

  /// \brief Identifier type for threads (i.e. with OpenMP).
  using ThreadId = int;

  /// \brief Identifier type for tasks (i.e. with OpenMP).
  using TaskId = int;
}

#endif // DATATYPES_H
