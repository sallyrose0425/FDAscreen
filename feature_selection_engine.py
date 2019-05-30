import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from deap import base
from deap import creator
from deap import tools
import warnings
warnings.simplefilter(action='ignore')

model = RandomForestClassifier()

features = np.random.random((100, 100))
labels = np.random.choice(2, size=(100, 1))

_, num_features = np.shape(features)


def mask_value(mask):
    """Returns the mean roc_auc_score of a random forest model trained with
    the indicated subset of features."""
    features_tmp = features[:, mask == 1]
    labels_tmp = labels
    np.random.seed(42)
    scores = []
    mini_batches_generator = StratifiedKFold(n_splits=4, random_state=42, shuffle=True)
    try:
        for training_index, validation_index in mini_batches_generator.split(features_tmp, labels_tmp):
            training_features = features_tmp[training_index]
            training_labels = np.ravel(labels_tmp[training_index])
            validation_features = features_tmp[validation_index]
            validation_labels = np.ravel(labels_tmp[validation_index])
            model.fit(training_features, training_labels)
            predictions = model.predict_proba(validation_features)[:, 1]
            scores.append(roc_auc_score(validation_labels, predictions))
            return np.mean(scores)
    except ValueError:
        return 0


def mask_opt_function(mask):
    """The function being maximized by genetic algorithm.
     It attempts to balance the output of mask_value with the feature dimension."""
    feature_dim = np.sum(mask)
    model_auc = mask_value(mask)
    return model_auc ** 2 / (1 + feature_dim),


def genetic_algorithm(num_gens):
    np.random.seed(42)
    print_freq = 10
    target_ratio = 0.1
    pop_size = 50
    INDPB = 0.05
    TOURNSIZE = 3
    CXPB = 0.5
    MUTPB = 0.2

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", np.ndarray, typecode='b', fitness=creator.FitnessMax)
    toolbox = base.Toolbox()
    toolbox.register("attr_bool", np.random.choice, 2, p=[1 - target_ratio, target_ratio])
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, num_features)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", mask_opt_function)
    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=INDPB)
    toolbox.register("select", tools.selTournament, tournsize=TOURNSIZE)

    max_record = []
    pop = toolbox.population(n=pop_size)
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    gen = 0
    while gen < num_gens:
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if np.random.random() < CXPB:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        for mutant in offspring:
            if np.random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        pop[:] = offspring
        scores = [mask_value(mask) for mask in pop]
        if gen % print_freq == 0:
            best_mask = pop[np.argmax(scores)]
            max_score = np.round(np.max(scores), 3)
            print('Best ROC-AUC: {},'
                  ' Feature Dimension: {}'.format(max_score, np.sum(best_mask)))
            max_record.append([max_score, np.sum(best_mask)])
        gen += 1
    return max_record, pop[np.argmax(scores)]


max_record, feature_mask = genetic_algorithm(1000)
