#include <TestUtils/base.h>

#include <Core/lm_vector.h>

using namespace Core;

TEST(LMVector, Works)
{
  const auto empty0 = lm_vector<double>(0);
  EXPECT_EQ(empty0.l_max, 0u);
  const auto empty3 = lm_vector<double>(3);
  EXPECT_EQ(empty3.l_max, 3u);

  const auto non_empty0 = lm_vector<double>(std::vector<double>{1.0});
  EXPECT_EQ(non_empty0.l_max, 0u);
  EXPECT_DOUBLE_EQ(non_empty0.at(0, 0), 1.0);
  const auto non_empty2 = lm_vector<double>(std::vector<double>{
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
  });
  EXPECT_EQ(non_empty2.l_max, 2u);
  EXPECT_DOUBLE_EQ(non_empty2.at(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(1, -1), 1.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(1, 1), 3.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(2, -2), 4.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(2, -1), 5.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(2, 0), 6.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(2, 1), 7.0);
  EXPECT_DOUBLE_EQ(non_empty2.at(2, 2), 8.0);
}

TEST(LMVector, Throws)
{
  const auto data = std::vector<double>{0.0, 1.0, 2.0};
  auto ctor = [&data]() -> void { auto a = lm_vector<double>(data); };
  EXPECT_THROW(ctor(), std::invalid_argument);

  const auto non_empty2 = lm_vector<double>(std::vector<double>{
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
  });
  EXPECT_THROW(non_empty2.at(1, -2), std::invalid_argument);
  EXPECT_THROW(non_empty2.at(0, 1), std::invalid_argument);
  EXPECT_THROW(non_empty2.at(2, 5), std::invalid_argument);
  EXPECT_THROW(non_empty2.at(3, 1), std::invalid_argument);
}