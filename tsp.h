/*************************************************************************
Title: TSP.hpp
Description: TSP class specification file for our Christofides implementation
Authors: Sean Hinds, Ryan Hong, Jeff Herlitz
Date: 08/16/17
*************************************************************************/
#ifndef TSP_H
#define TSP_H
//---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stack>
#include <string>
#include <stdio.h>
#include <vector>
#include <limits>
//---------------------------------------------------------------------------
//1=int, 2=float, 3=double, 4=long double
#define NUM_TYPE 3

#if NUM_TYPE==1
#define TSPN_NUM_TYPE int
#elif NUM_TYPE==2
#define TSPN_NUM_TYPE float
#elif NUM_TYPE==3
#define TSPN_NUM_TYPE double
#elif NUM_TYPE==4
#define TSPN_NUM_TYPE long double
#endif
//---------------------------------------------------------------------------
struct TSPN_Point{
    TSPN_NUM_TYPE x;
    TSPN_NUM_TYPE y;
};
//---------------------------------------------------------------------------
class TSPN
{
 protected:
  // List of odd nodes
  std::vector<int> odds;

  //Adjacency list
  std::vector<int>* adjlist;

  // Number of points
  int n;

  //Shortest path length
  TSPN_NUM_TYPE pathLength;

  //euler circuit
  std::vector<int> circuit;

  // n x n, pairwise distances between points
  TSPN_NUM_TYPE **graph;

  // Point list
  std::vector<TSPN_Point> points;

  // -
  void findOdds();

  //Find perfect matching
  void perfectMatching();

  //Find Euler tour
  void euler_tour(int start, std::vector<int> &path);
  
  //Find Hamiltonian path
  void make_hamiltonian(std::vector<int> &path, TSPN_NUM_TYPE &pathCost);

  TSPN_NUM_TYPE get_distance(struct TSPN_Point c1, struct TSPN_Point c2);
  void findMST();
  int getMinIndex(TSPN_NUM_TYPE key[], bool mst[]);
  void fillMatrix();
  TSPN_NUM_TYPE findBestPath(int start);

 public:

  TSPN(std::vector<TSPN_Point> aPointList);
  ~TSPN();

  void Solve();

  void printResult();
  void printPath();
  void printAdjList();

  std::vector<int> getResult();
};
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
