#pragma once

/**
 *@brief Enum class for boundary conditions
 */

enum BoundaryCondition { Outflow, Reflecting };

/**
 *@brief Struct for domain boundary conditions
 */

struct DomainBoundaryConditions {
  BoundaryCondition left, right;
  BoundaryCondition top, bottom;
  BoundaryCondition front, back;
};