
## mgatk python dev notes
**Caleb Lareau**

### DAG

To visualize the DAG associated with the `snakemake` file, do the following--

```
1) Modify the Snakefile.Scatter to 
2) Run a basic mgatk command (e.g. mgatk call -i tests/humanbam -o tests/outbam)
```

This will generate some plaintext that looks something like this:

```
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "all", color = "0.40 0.6 0.85", style="rounded"];
	1[label = "make_depth_table", color = "0.13 0.6 0.85", style="rounded"];
	2[label = "call_variants", color = "0.53 0.6 0.85", style="rounded"];
	3[label = "get_depth_bq_baq", color = "0.27 0.6 0.85", style="rounded"];
	4[label = "filter_sort", color = "0.00 0.6 0.85", style="rounded"];
	1 -> 0
	2 -> 0
	3 -> 1
	4 -> 2
	4 -> 3
}     
```

So dump those contents into a text file, say, `dag.txt`. Finally, run--

```
cat dag.txt | dot -Tsvg > dag.svg
```
to visualize the dag. 
<br><br>