import java.util.*;
import java.util.ArrayList;
import java.util.List;

public class EvolutionAlgorithm {
    private final int parentsNumber;
    private final int childsPerParentNumber;
    private final int epochsNumber;
    List<double[]> population;
    List<Double> populationScores;
    private final FitnessFunction fitnessFunction;
    private final double mutationRate;
    private final double sigma; //середньоквадратичне відхилення
    private long evolutionDuration;
    private final Random random;
    private double bestPopulationScore;

    public EvolutionAlgorithm(int parentsNumber, int childsPerParentNumber, int epochsNumber, double mutationRate, double sigma,
                              FitnessFunction fitnessFunction, double[] bounds) {
        this.random = new Random();
        this.parentsNumber = parentsNumber;
        this.childsPerParentNumber = childsPerParentNumber;
        this.epochsNumber = epochsNumber;
        this.fitnessFunction = fitnessFunction;
        this.mutationRate = mutationRate;
        this.sigma = sigma;
        population = new ArrayList<>();
        populationScores = new ArrayList<>();
        createStartPopulation(bounds);
    }

    public void startEvolution() {
        long startTime = System.currentTimeMillis();
        for (int epoch = 0; epoch < epochsNumber; ++epoch) {
            calculatePopulationScores();
            addChilds();
            calculatePopulationScores();
            formNewPopulation();
            mutateAllIndivids();
        }

        evolutionDuration = System.currentTimeMillis() - startTime;
    }

    private void createStartPopulation(double[] bounds) {
        double mean = (bounds[0] + bounds[1]) / 2.0;
        double stdDev = (bounds[1] - bounds[0]) / 6;
        for (int i = 0; i < parentsNumber; ++i) {
            this.population.add(new double[]{random.nextGaussian()* stdDev + mean, random.nextGaussian()* stdDev + mean});
        }
    }

    private void addChilds() {
        List<double[]> childs = new ArrayList<>();

        //(lambda, mu)//
/*        for (int i = 0; i < parentsNumber; ++i) {
            double[] parent = population.get(random.nextInt(population.size()));
            for (int j = 0; j < childsPerParentNumber; ++j) {
                double[] child = Arrays.copyOf(parent, parent.length);
                child[0] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                child[1] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                childs.add(child);
            }
        }*/

        //(lambda/ro, mu)//
        for (int i = 0; i < parentsNumber/2; ++i){
            for (int j = 0; j < childsPerParentNumber; ++j){
                double[] parent1 = tournamentSelection();
                double[] parent2 = tournamentSelection();
                childs.add(new double[]{ (parent1[0] + parent2[0])/2, (parent1[1] + parent2[1])/2 });
            }
        }

        //(lambda, mu)//
        //population.clear();
        population.addAll(childs);
    }

    private void formNewPopulation() {
        List<double[]> newPopulation = new ArrayList<>(population);

        newPopulation.sort((o1, o2) -> Double.compare(fitnessFunction.calculate(o1[0], o1[1]), fitnessFunction.calculate(o2[0], o2[1])));

        List<double[]> selectedPopulation = newPopulation.subList(0, parentsNumber);
        bestPopulationScore = fitnessFunction.calculate(selectedPopulation.get(0)[0], selectedPopulation.get(0)[1]);

        population.clear();
        population.addAll(selectedPopulation);
    }

    private void mutateAllIndivids() {
        for (double[] doubles : population)
            if (random.nextDouble() < mutationRate) {
                doubles[0] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                doubles[1] += scaleToRange(random.nextGaussian(), -sigma, sigma);
            }
    }

    private double scaleToRange(double value, double min, double max) { return ((value + 1.0) / 2.0) * (max - min) + min; }

    private void calculatePopulationScores() {
        if (!populationScores.isEmpty())
            populationScores.clear();

        for (double[] individ : population)
            populationScores.add(fitnessFunction.calculate(individ[0], individ[1]));
    }

    private double[] tournamentSelection() {
        double[] bestIndividual = null;
        double bestScore = Double.MAX_VALUE;

        for (int i = 0; i < 3; i++) {
            int randomIndex = random.nextInt(population.size());
            double score = populationScores.get(randomIndex);

            if (score < bestScore) {
                bestIndividual = population.get(randomIndex);
                bestScore = score;
            }
        }
        return bestIndividual;
    }

    public double getBestScore() { return bestPopulationScore; }
    public long getEvolutionDuration() { return evolutionDuration; }
}
