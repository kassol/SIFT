#include "match.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
	match m("D:\\01-155_50mic.tif", "D:\\01-156_50mic.tif");
	std::vector<SamePoint> result;
	m.domatch(result);
	system("pause");
	return 0;
}
