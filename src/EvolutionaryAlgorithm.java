import java.util.*;

public class EvolutionaryAlgorithm {
    private final int lambda;
    private final int mu;
    private final int epochsNum;
    List<double[]> population;
    List<Double> populationScores;
    private final FitnessFunction fitnessFunction;
    private final double mutationRate;
    private final double sigma; //середньоквадратичне відхилення
    private long evolutionDuration;
    private final Random random;
    private double bestPopulationScore;

    public EvolutionaryAlgorithm(int lambda, int mu, int epochsNum, double mutationRate, double sigma, FitnessFunction fitnessFunction, double[] bounds) {
        this.random = new Random();
        this.lambda = lambda;
        this.mu = mu;
        this.epochsNum = epochsNum;
        this.fitnessFunction = fitnessFunction;
        this.mutationRate = mutationRate;
        this.sigma = sigma;
        population = new ArrayList<>();
        populationScores = new ArrayList<>();
        createPopulation(bounds);
    }

    public void startEvolution() {
        long startTime = System.currentTimeMillis();
        for (int epoch = 0; epoch < epochsNum; ++epoch) {
            calculateScores();
            addChilds();
            calculateScores();
            selectNextParents();
            mutate();
        }

        evolutionDuration = System.currentTimeMillis() - startTime;
    }

    private void createPopulation(double[] bounds) {
        for (int i = 0; i < lambda; ++i) 
            this.population.add(new double[]{random.nextGaussian()* (bounds[1] - bounds[0]) / 6 + (bounds[0] + bounds[1]) / 2.0, random.nextGaussian() * (bounds[1] - bounds[0]) / 6 + (bounds[0] + bounds[1]) / 2.0});
    }

    private void addChilds() {
        strategyLambdaPhoMu();
        //strategyLambdaPlusMu();
        //strategyLambdaMu();
    }

    private void strategyLambdaPlusMu() {
        List<double[]> childs = new ArrayList<>();

        for (int i = 0; i < lambda; ++i) {
            double[] parent = population.get(random.nextInt(population.size()));
            for (int j = 0; j < mu; ++j) {
                double[] child = Arrays.copyOf(parent, parent.length);
                child[0] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                child[1] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                childs.add(child);
            }
        }

        population.addAll(childs);
    }

    private void strategyLambdaPhoMu () {
        List<double[]> childs = new ArrayList<>();

        for (int i = 0; i < lambda/2; ++i){
            for (int j = 0; j < mu; ++j){
                double[] parent1 = tournamentSelection();
                double[] parent2 = tournamentSelection();
                childs.add(new double[]{ (parent1[0] + parent2[0])/2, (parent1[1] + parent2[1])/2 });
            }
        }

        population.addAll(childs);
    }

    private void strategyLambdaMu () {
        List<double[]> childs = new ArrayList<>();

        for (int i = 0; i < lambda; ++i) {
            double[] parent = population.get(random.nextInt(population.size()));
            for (int j = 0; j < mu; ++j) {
                double[] child = Arrays.copyOf(parent, parent.length);
                child[0] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                child[1] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                childs.add(child);
            }
        }

        population.clear();
        population.addAll(childs);
    }

    private void selectNextParents() {
        List<double[]> nextPopulation = new ArrayList<>(population);

        nextPopulation.sort(Comparator.comparingDouble(individual -> fitnessFunction.calculate(individual[0], individual[1])));

        nextPopulation = nextPopulation.subList(0, lambda);
        bestPopulationScore = fitnessFunction.calculate(nextPopulation.get(0)[0], nextPopulation.get(0)[1]);

        population.clear();
        population.addAll(nextPopulation);
    }

    private void mutate() {
        for (double[] doubles : population)
            if (random.nextDouble() < mutationRate) {
                doubles[0] += scaleToRange(random.nextGaussian(), -sigma, sigma);
                doubles[1] += scaleToRange(random.nextGaussian(), -sigma, sigma);
            }
    }

    private double scaleToRange(double value, double min, double max) {
        return ((value + 1.0) / 2.0) * (max - min) + min;
    }

    private void calculateScores() {
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
