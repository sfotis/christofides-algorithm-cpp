/*************************************************************************
Title: TSP.hpp
Description: TSP class specification file for our Christofides implementation
Authors: Sean Hinds, Ryan Hong, Jeff Herlitz
Date: 08/16/17
/*************************************************************************
Modified & enhanced by Fotios Sioutis "sfotis@gmail.com"
Date: 06/11/19
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
template <typename T> class TSP
{
 public:
  struct TSP_Point{
    T x;
    T y;
  };

 protected:
  // List of odd nodes
  std::vector<int> odds;

  //Adjacency list
  std::vector<int>* adjlist;

  // Number of points
  int n;

  //Shortest path length
  T pathLength;

  //euler circuit
  std::vector<int> circuit;

  // n x n, pairwise distances between points
  T **graph;

  // Point list
  std::vector<TSP_Point> points;

  // -
  void findOdds();

  //Find perfect matching
  void perfectMatching();

  //Find Euler tour
  void euler_tour(int start, std::vector<int> &path);
  
  //Find Hamiltonian path
  void make_hamiltonian(std::vector<int> &path, T &pathCost);

  T get_distance(struct TSP_Point c1, struct TSP_Point c2);

  //T get_distance(struct TSP_Point c1, struct TSP_Point c2);

  void findMST();
  int getMinIndex(T key[], bool mst[]);
  void fillMatrix();
  T findBestPath(int start);

 public:
  TSP(std::vector<TSP_Point> aPointList);
  ~TSP();

  void Solve();

  void printResult();
  void printPath();
  void printAdjList();

  std::vector<int> getResult();
};
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
