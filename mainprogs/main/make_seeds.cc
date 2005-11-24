/*
 *     This is rannyu as modified by A. Sokal 9/26/85.
 *        It is linear congruential with modulus m = 2**48, increment c = 1,
 *        and multiplier a = (2**36)*m1 + (2**24)*m2 + (2**12)*m3 + m4. 
 *        The multiplier is stored in common (see subroutine setrn)
 *      and is set to a = 31167285 (recommended by Knuth, vol. 2,
 *      2nd ed., p. 102).
 *
 *     Multiplier is 31167285 = (2**24) + 3513*(2**12) + 821.
 *        Recommended by Knuth, vol. 2, 2nd ed., p. 102.
 *     (Generator is linear congruential with odd increment
 *        and maximal period, so seed is unrestricted: it can be
 *        either even or odd.)
 */

static int m[4] = {0, 1, 3513, 821};
static int l[4] = {0, 0, 0, 1};

float rannyu()
{
  float twom12 = 1/4096.0;
  int i[4];

  i[0] = l[0]*m[3] + l[1]*m[2] + l[2]*m[1] + l[3]*m[0];
  i[1] = l[1]*m[3] + l[2]*m[2] + l[3]*m[1];
  i[2] = l[2]*m[3] + l[3]*m[2];
  i[3] = l[3]*m[3] + 1;
  l[3] = i[3] & 4095;
  i[2] = i[2] + (i[3] >> 12);
  l[2] = i[2] & 4095;
  i[1] = i[1] + (i[2] >> 12);
  l[1] = i[1] & 4095;
  l[0] = (i[0] + (i[1] >> 12)) >> 12;
  return twom12*((float)l[0] + twom12*((float)l[1] + twom12*((float)l[2] + twom12*((float)l[3]))));
}

/*! Seed has been set by default - this allows one to override it */
void setrn(int iseed[4])
{
  int i;
  for(i=0; i < 4; ++i)
    l[i] = iseed[i];
}

/*! Recover the seed */
void savern(int iseed[4])
{
  int i;
  for(i=0; i < 4; ++i)
    iseed[i] = l[i];
}


#include <cstdio>
#include <iostream>

int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage:  " << argv[0] << ": <num throw away> <num seeds>\n";
    exit(1);
  }

  int ndisc = atoi(argv[1]);
  int num = atoi(argv[2]);

  for(int n=0; n < ndisc; ++n)
    rannyu();

  for(int n=0; n < num; ++n)
  {
    float f = rannyu();

    int iseed[4];
    savern(iseed);
    printf("%d %d %d %d\n", iseed[1], iseed[2], iseed[3], iseed[0] & 2047);
  }
}
