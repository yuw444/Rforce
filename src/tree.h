/*
@author : Yu Wang
*/

#ifndef TREE_H
#define TREE_H
#include "utils.h"
#include "split.h"
#include <omp.h>

// global variables
// flag for dynamic risk time estimation at each bootstrap iteration; 1 is on; 0 is off
// verbose for debugging; 1 is on; 0 is off
extern int verbose;

// struct to store the structure of one node
typedef struct DecisionTreeNode
{
    long treeId;
    long nodeId;

    unsigned int flag; // indicate if it's a terminal node, 1 is terminal node
    long splitIndex;   // only based on designMatrix only
    double splitValue; // only based on designMatrix only
    double *splitStat; // could use anxiliaryFeatures;
    size_t *sizeLR;    // the size of daughter nodes;
    double *output;    // could use auxiliaryFeatures;
    size_t lenOutput; // the length of output

    struct DecisionTreeNode *leftChild;
    struct DecisionTreeNode *rightChild;

} DecisionTreeNode;

// create a new node
DecisionTreeNode *EmptyNode(long treeId, long *nodeId);

// grow a tree
void GrowTree(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    DecisionTreeNode *tree,
    double **designMatrixY,     // only a copy of pointer to rows, not entire data matrix
    double **auxiliaryFeatures, // only a copy of pointer to rows, not entire data matrix
    size_t nrows,               // nrow for design matrix
    size_t ncolsDesign,         // ncol for design matrix, including Y
    size_t ncolsAuxiliary,      // ncol for auxiliaryFeature matrix
    double *unitsOfCPIU,        // length of each CPIU
    size_t nUnits,              // number of CPIUs
    size_t lenOutput,           // length of output
    long treeId,               // current tree id
    long *nodeId,               // current node id
    size_t depth,               // current depth
    size_t maxDepth,            // max depth
    size_t minNodeSize,         // min size of leaf node
    double minGain,             // min gain of split
    size_t mtry,                // mtry variables when splitting
    size_t nsplits,             // maximum attempts of each variable when splitting
    unsigned int seed);         // seed for random number generator

// single tree functions; A wrapper of GrowTree();
DecisionTreeNode *Tree(
    fpSplitFunction splitFunction,
    fpLeafOutput leafOutputFunction,
    double **designMatrixY,
    double **auxiliaryFeatures,
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    double *unitsOfCPIUs,
    size_t nUnits,
    size_t lenOutput,
    long treeId,
    size_t maxDepth,
    size_t minNodeSize,
    double minGain,
    size_t mtry,
    size_t nsplits,
    unsigned int seed
);

// copy a tree
DecisionTreeNode *TreeCopy(DecisionTreeNode *root);

// recursive function to 
// 1. partition the data
// 2. prune the tree structure
// 3. calculate the quasilikelihood at the terminal nodes
void TreeOOBInternal(
    DecisionTreeNode *rootCopy,
    double **designMatrixY,     // the out-of-bag design matrix
    double **auxiliaryFeatures, // the out-of-bag auxiliary features, pseudo-risk estimation should be done in forest step before feeding into this function, similar with GrowTree()
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    size_t lenOutput);

// wrapper function for TreeOOBInternal to reduce the code complexity by saving the deference for pointer
DecisionTreeNode *TreeOOB(
    DecisionTreeNode *root,
    double **designMatrixY,     // the out-of-bag design matrix
    double **auxiliaryFeatures, // the out-of-bag auxiliary features, pseudo-risk estimation should be done in forest step before feeding into this function, similar with GrowTree()
    size_t nrows,
    size_t ncolsDesign,
    size_t ncolsAuxiliary,
    size_t lenOutput);

double TreeLikelihoodSum(DecisionTreeNode *root);

// free memory for a tree
void FreeTree(
    DecisionTreeNode *root);

// predict for a single row of testing data
double *TreePredict(
    DecisionTreeNode *root,
    double *CPIU); // one row of testing data 
    
// print tree in dot format
void PrintTreeDot(
    FILE *file,
    DecisionTreeNode *root,
    unsigned int level);

// A wrapper of PrintTreeDot() and write to file
void WriteTreeDotFile(
    DecisionTreeNode *root,
    char *filename,
    unsigned int level);

// gather vimp vector for a tree
void VimpTree(
    DecisionTreeNode *root,
    double **vimpStat,
    double **vimpFreq);

// print tree for debugging
void PrintTree(DecisionTreeNode *root, unsigned int level);

// save and load tree
void SaveTree(
    DecisionTreeNode *root,
    FILE *file);

DecisionTreeNode *LoadTree(FILE *file);
// print all tree nodes elements

typedef struct NodeElements
{
    long treeId;
    long nodeId;

    unsigned int flag;
    long splitIndex;
    double splitValue;
    double *splitStat;
    size_t *sizeLR;
    double *output;
    size_t lenOutput;

} NodeElements;

void PrintAllNodesElementsRecur(
    DecisionTreeNode *root,
    FILE *file,
    NodeElements **nodeElements,
    size_t pathDepth,
    size_t *nthPath);

void PrintTreeNodes(
    DecisionTreeNode *root,
    FILE *file);

#endif