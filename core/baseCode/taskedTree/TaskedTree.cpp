/*
 * file: TaskedTree.cpp
 * description: TaskedTree processing package.
 * author: Gueunet Charles
 * date: Dec 2016
 */

#include "TaskedTree.h"

// ------------
// CONSTRUCTOR
// ------------
// {

TaskedTree::TaskedTree() : ContourTree(new Params, nullptr, new Scalars)
{
}

TaskedTree::~TaskedTree()
{
}
