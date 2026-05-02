/* addtest_basic.C
 * Demonstrates usin the function imported from the DLL, the inelegant way
 */
#include <stdlib.h>
#include <stdio.h>

__declspec(dllimport) int __cdecl Add(int a, int b);

int main(int argc, char** argv) {
   printf("%d\n", Add(6, 23));

   return EXIT_SUCCESS;
}
