# Visualize a Single Tree from an Rforce Forest

This function extracts a single decision tree from an `Rforce` forest,
generates its Graphviz DOT representation via the underlying C backend,
and renders it using DiagrammeR. The DOT file is stored temporarily and
automatically removed after the R session.

## Usage

``` r
printTree(forest, treeIndex)
```

## Arguments

- forest:

  An object of class `Rforce` containing the trained forest and a valid
  external C pointer in `` `_external_forest_C_Ptr` ``.

- treeIndex:

  An integer specifying the index of the tree to visualize. Must be
  between `1` and `forest\$numTrees`.

## Value

A DiagrammeR graph object representing the visualized decision tree.
Additionally, a message is printed showing the location of the temporary
DOT file.

## Details

Internally, this function:

1.  Validates input.

2.  Creates a temporary `.dot` file.

3.  Calls the C function `R_PrintTree` via `.Call` to write the DOT
    graph representation of the tree.

4.  Displays the tree structure using
    [`DiagrammeR::grViz()`](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html).

The DOT file persists only for the current R session and is stored in a
temporary directory. It is automatically cleaned up when the session
ends.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming `rf_model` is an Rforce object with at least one tree:
printTree(rf_model, treeIndex = 1)
} # }
```
