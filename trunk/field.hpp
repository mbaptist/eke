


template <class Grid,class Value>
class Field
{
private:
  Grid grid_;
  Value value_;
public:
  Grid & grid(){return grid_;}; 
  Value & value(){return value_;};
  const Grid & grid() const {return grid_;}; 
  const Value & value() const {return value_;};
public:
  Field(Grid _grid_,Value _value_):
    grid_(_grid_),
    value_(_value_)
  {
  };
  ~Field(){};
private:
  Field();
  
};

