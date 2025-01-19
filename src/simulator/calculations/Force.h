#include "../particle/container/DirectSumContainer.h"
#include "../particle/container/LinkedCellContainer.h"
#include "Calculation.h"

/**
 * @enum ForceType
 * @brief Enum class for the force calculation type.
 */
enum ForceType { LENNARD_JONES, GRAVITATIONAL };

/***
 * @struct Force
 * @brief Struct, that provides functions for force calculation.
 **/
struct Force : AbstractPolicy {
  static void run(LinkedCellContainer &particles, ForceType type,
                  OPTIONS OPTION);

private:
  /***
   * @brief Lennard-Jones force calculation.
   * @param particles DirectSumContainer reference.
   **/
  static void lennard_jones(LinkedCellContainer &particles, OPTIONS OPTION);
  /**
   * @brief Gravitational force calculation.
   * @param particles DirectSumContainer reference.
   */
  static void gravitational(LinkedCellContainer &particles, OPTIONS OPTION);
};