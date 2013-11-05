#include "match.h"
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
	match m("D:\\left.png", "D:\\right.png");
	std::vector<SamePoint> result;
	m.domatch(result);
	system("pause");
	return 0;
}
