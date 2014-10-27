
#include <vector>
#include <algorithm>

struct S
{
int field1;
double field2;
};

struct SortByField2
{
bool operator () (const S & lhs , const S & rhs) const
{
return lhs.field2 < rhs.field2;
}
};

int main()
{
std::vector<S> v;

// add to vector

std::sort(v.begin(),v.end(),SortByField2());

return 0;
}
