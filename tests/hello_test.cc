#include <gtest/gtest.h>
#include "ParticleContainer.h"

TEST(HelloTest, BasicAssertions) {
  ParticleContainer p{};
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
}