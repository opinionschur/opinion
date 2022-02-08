# Maximizing Opinion in Social Networks via Leader Selection

Julia code for the article "Maximizing Opinion in Social Networks via Leader Selection"


[edgecore.jl](./edgecore.jl) and [alg1.jl](./alg1.jl) contain Greedy and Fast and other baselines mentioned in the article, which are function exactgreedy() and function greedy() and so on.


[edgemain.jl](./edgemain.jl) is the main body of the algorithm and calls function exactgreedy() and function greedy() and other baselines to compare their effectiveness and efficiency


[soc-dolphins.mtx](./data/soc-dolphins.mtx) is an example of the data
