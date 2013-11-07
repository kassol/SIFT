#include "match.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
	match m("D:\\1.tif", "D:\\2.tif");
	std::vector<SamePoint> result;
	m.domatch(result);
	system("pause");
	return 0;
}
