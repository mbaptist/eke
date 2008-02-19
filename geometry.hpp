
#include <blitz/array.h>

class GeometryObject
{
};


class Sphere: public GeometryObject
{
  
};


class Cube: public GeometryObject
{
};


class Geometry: public blitz::Array<GeometryObject,1>
{
};

