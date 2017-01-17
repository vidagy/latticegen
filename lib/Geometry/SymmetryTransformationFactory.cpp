#include <algorithm>
#include <set>
#include <iostream>
#include "SymmetryTransformationFactory.h"

namespace Geometry
{
  namespace
  {
    static const double pi = 3.14159265358979323846;

    static const Vector3D x = {1.0, 0.0, 0.0};
    static const Vector3D y = {0.0, 1.0, 0.0};
    static const Vector3D z = {0.0, 0.0, 1.0};

    static const Vector3D p30 = {sqrt(3.0)/2.0, 0.5, 0.0};
    static const Vector3D m30 = {sqrt(3.0)/2.0, -0.5, 0.0};
    static const Vector3D p60 = {0.5, sqrt(3.0)/2.0, 0.0};
    static const Vector3D m60 = {0.5, -sqrt(3.0)/2.0, 0.0};

    static const Vector3D p1p10 = {1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0};
    static const Vector3D p1m10 = {1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0.0};
    static const Vector3D p10p1 = {1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0)};
    static const Vector3D p10m1 = {1.0 / sqrt(2.0), 0.0, -1.0 / sqrt(2.0)};
    static const Vector3D Np1p1 = {0.0, 1.0 / sqrt(2.0), 1.0 / sqrt(2.0)};
    static const Vector3D Np1m1 = {0.0, 1.0 / sqrt(2.0), -1.0 / sqrt(2.0)};

    static const Vector3D p1p1p1 = 1.0 / sqrt(3.0) * Vector3D{ 1.0, 1.0, 1.0};
    static const Vector3D p1m1m1 = 1.0 / sqrt(3.0) * Vector3D{ 1.0,-1.0,-1.0};
    static const Vector3D m1p1m1 = 1.0 / sqrt(3.0) * Vector3D{-1.0, 1.0,-1.0};
    static const Vector3D m1m1p1 = 1.0 / sqrt(3.0) * Vector3D{-1.0,-1.0, 1.0};
  }


  Transformation SymmetryTransformationFactory::get(
    const CrystallographicPointGroup::SymmetryElement symmetry_element)
  {
    switch (symmetry_element)
    {
      case CrystallographicPointGroup::Identity:
        return Identity();
      case CrystallographicPointGroup::Inversion:
        return Inversion();

      case CrystallographicPointGroup::Rotation2100:
        return Rotation(pi * x);
      case CrystallographicPointGroup::Rotation2010:
        return Rotation(pi * y);
      case CrystallographicPointGroup::Rotation2001:
        return Rotation(pi * z);
      case CrystallographicPointGroup::Rotation2p300:
        return Rotation(pi * p30);
      case CrystallographicPointGroup::Rotation2m300:
        return Rotation(pi * m30);
      case CrystallographicPointGroup::Rotation2p600:
        return Rotation(pi * p60);
      case CrystallographicPointGroup::Rotation2m600:
        return Rotation(pi * m60);
      case CrystallographicPointGroup::Rotationp3001:
        return Rotation(2.0 / 3.0 * pi * z);
      case CrystallographicPointGroup::Rotationm3001:
        return Rotation(-2.0 / 3.0 * pi * z);
      case CrystallographicPointGroup::Rotationp4001:
        return Rotation(pi / 2.0 * z);
      case CrystallographicPointGroup::Rotationm4001:
        return Rotation(-pi / 2.0 * z);
      case CrystallographicPointGroup::Rotationp4010:
        return Rotation(pi / 2.0 * y);
      case CrystallographicPointGroup::Rotationm4010:
        return Rotation(-pi / 2.0 * y);
      case CrystallographicPointGroup::Rotationp4100:
        return Rotation(pi / 2.0 * x);
      case CrystallographicPointGroup::Rotationm4100:
        return Rotation(-pi / 2.0 * x);
      case CrystallographicPointGroup::Rotationp6001:
        return Rotation(pi / 3.0 * z);
      case CrystallographicPointGroup::Rotationm6001:
        return Rotation(-pi / 3.0 * z);
      case CrystallographicPointGroup::Rotation2p1p10:
        return Rotation(pi * p1p10);
      case CrystallographicPointGroup::Rotation2p1m10:
        return Rotation(pi * p1m10);
      case CrystallographicPointGroup::Rotation2p10p1:
        return Rotation(pi * p10p1);
      case CrystallographicPointGroup::Rotation2p10m1:
        return Rotation(pi * p10m1);
      case CrystallographicPointGroup::Rotation20p1p1:
        return Rotation(pi * Np1p1);
      case CrystallographicPointGroup::Rotation20p1m1:
        return Rotation(pi * Np1m1);
      case CrystallographicPointGroup::Rotationp3p1p1p1:
        return Rotation(2.0 / 3.0 * pi * p1p1p1);
      case CrystallographicPointGroup::Rotationm3p1p1p1:
        return Rotation(-2.0 / 3.0 * pi * p1p1p1);
      case CrystallographicPointGroup::Rotationp3p1m1m1:
        return Rotation(2.0 / 3.0 * pi * p1m1m1);
      case CrystallographicPointGroup::Rotationm3p1m1m1:
        return Rotation(-2.0 / 3.0 * pi * p1m1m1);
      case CrystallographicPointGroup::Rotationp3m1p1m1:
        return Rotation(2.0 / 3.0 * pi * m1p1m1);
      case CrystallographicPointGroup::Rotationm3m1p1m1:
        return Rotation(-2.0 / 3.0 * pi * m1p1m1);
      case CrystallographicPointGroup::Rotationp3m1m1p1:
        return Rotation(2.0 / 3.0 * pi * m1m1p1);
      case CrystallographicPointGroup::Rotationm3m1m1p1:
        return Rotation(-2.0 / 3.0 * pi * m1m1p1);

      case CrystallographicPointGroup::Reflection100:
        return Reflection(x);
      case CrystallographicPointGroup::Reflection010:
        return Reflection(y);
      case CrystallographicPointGroup::Reflection001:
        return Reflection(z);
      case CrystallographicPointGroup::Reflectionp1p10:
        return Reflection(p1p10);
      case CrystallographicPointGroup::Reflectionp1m10:
        return Reflection(p1m10);
      case CrystallographicPointGroup::Reflectionp10p1:
        return Reflection(p10p1);
      case CrystallographicPointGroup::Reflectionp10m1:
        return Reflection(p10m1);
      case CrystallographicPointGroup::Reflection0p1p1:
        return Reflection(Np1p1);
      case CrystallographicPointGroup::Reflection0p1m1:
        return Reflection(Np1m1);
      case CrystallographicPointGroup::Reflectionp300:
        return Reflection(p30);
      case CrystallographicPointGroup::Reflectionm300:
        return Reflection(m30);
      case CrystallographicPointGroup::Reflectionp600:
        return Reflection(p60);
      case CrystallographicPointGroup::Reflectionm600:
        return Reflection(m60);

      case CrystallographicPointGroup::ImproperRotationp3001:
        return ImproperRotation(2.0 / 3.0 * pi * z);
      case CrystallographicPointGroup::ImproperRotationm3001:
        return ImproperRotation(-2.0 / 3.0 * pi * z);
      case CrystallographicPointGroup::ImproperRotationp4001:
        return ImproperRotation(pi / 2.0 * z);
      case CrystallographicPointGroup::ImproperRotationm4001:
        return ImproperRotation(-pi / 2.0 * z);
      case CrystallographicPointGroup::ImproperRotationp4010:
        return ImproperRotation(pi / 2.0 * y);
      case CrystallographicPointGroup::ImproperRotationm4010:
        return ImproperRotation(-pi / 2.0 * y);
      case CrystallographicPointGroup::ImproperRotationp4100:
        return ImproperRotation(pi / 2.0 * x);
      case CrystallographicPointGroup::ImproperRotationm4100:
        return ImproperRotation(-pi / 2.0 * x);
      case CrystallographicPointGroup::ImproperRotationp6001:
        return ImproperRotation(pi / 3.0 * z);
      case CrystallographicPointGroup::ImproperRotationm6001:
        return ImproperRotation(-pi / 3.0 * z);
      case CrystallographicPointGroup::ImproperRotationp6p1p1p1:
        return ImproperRotation(pi / 3.0 * p1p1p1);
      case CrystallographicPointGroup::ImproperRotationm6p1p1p1:
        return ImproperRotation(-pi / 3.0 * p1p1p1);
      case CrystallographicPointGroup::ImproperRotationp6p1m1m1:
        return ImproperRotation(pi / 3.0 * p1m1m1);
      case CrystallographicPointGroup::ImproperRotationm6p1m1m1:
        return ImproperRotation(-pi / 3.0 * p1m1m1);
      case CrystallographicPointGroup::ImproperRotationp6m1p1m1:
        return ImproperRotation(pi / 3.0 * m1p1m1);
      case CrystallographicPointGroup::ImproperRotationm6m1p1m1:
        return ImproperRotation(-pi / 3.0 * m1p1m1);
      case CrystallographicPointGroup::ImproperRotationp6m1m1p1:
        return ImproperRotation(pi / 3.0 * m1m1p1);
      case CrystallographicPointGroup::ImproperRotationm6m1m1p1:
        return ImproperRotation(-pi / 3.0 * m1m1p1);
      default:
        throw std::invalid_argument("symmetry_element not recognized in SymmetryTransformationFactory");
    }
  }

  SymmetryTransformationFactory::Transformations SymmetryTransformationFactory::get(
    const CrystallographicPointGroup::Elements elements)
  {
    Transformations transformations;
    transformations.reserve(elements.size());

    std::transform(elements.begin(), elements.end(), std::back_inserter(transformations),
                   [](const CrystallographicPointGroup::SymmetryElement symmetryElement) {
                     return SymmetryTransformationFactory::get(symmetryElement).transformation_matrix;;
                   }
    );
    return transformations;
  }

  SymmetryTransformationFactory::Transformations SymmetryTransformationFactory::generate(
    const CrystallographicPointGroup::Elements& generators)
  {
    Transformations new_elements = get(generators);
    Transformations elements{Identity()};
    Transformations buffer;
    buffer.reserve(48);

    while (new_elements.size())
    {
      Matrix3D new_element = new_elements.back();
      for(auto it = elements.begin(); it != elements.end(); ++it)
      {
        Matrix3D new_transformation = *it * new_element;
        auto it2 = std::find_if(
          elements.begin(),
          elements.end(),
          [&new_transformation](const Matrix3D &trafo)
          {
            return trafo == new_transformation;
          });
        if (it2 == elements.end()) {
          buffer.push_back(new_transformation);
        }
      }
      if (buffer.empty())
        new_elements.pop_back();
      else
      {
        std::copy(buffer.begin(), buffer.end(), std::back_inserter(new_elements));
        std::copy(buffer.begin(), buffer.end(), std::back_inserter(elements));
        buffer.clear();
      }
    }

    return Transformations{elements.begin(), elements.end()};
  }
}