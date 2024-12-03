#pragma once

enum class BoundaryCondition {
    Outflow,
    Reflecting
};

struct DomainBoundaryConditions {
    BoundaryCondition left, right;
    BoundaryCondition top, bottom;
    BoundaryCondition front, back;
};