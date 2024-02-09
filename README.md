# Improvement of ant colony algorithm based on block and search-Traveling salesman problem

We have improved the modified ant colony algorithm based on naive pheromone diffusion, and the effect of pheromone diffusion on the optimization of results has been improved through the block method, and the convergence rate can be controlled by adjusting the size of the block, the improvement code was based on C++

After applying our algorithm, we get a better result than the Naive Diffusion ACO algorithm, showing as the following table:

![alt text](reslut.png)

`ACO.ipynb` is the naive ACO algorithm, `NaiveDiffusionACO.ipynb` is the naive ACO algorithm based on the colony diffusion, and `BSACO.cpp` is the code we implemented.

`dantzig42.tsp` is the data we use for the ACO algorithm.

By using our algorithm, the partition result is shown in `partition.png` file.