#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.


from collections import defaultdict
from itertools import chain, cycle
from operator import itemgetter
from random import sample, choice

def selNSGA2(individuals, k):
    k = min(len(individuals), k)
    unique_fits = defaultdict(list)
    for ind in individuals:
        unique_fits[ind.fitness.wvalues].append(ind)
    fitnesses = unique_fits.keys()

    pareto_fronts, ranks = sortFastND(fitnesses, k)
    crowding_dist = computeCrowdingDist(pareto_fronts)
    chosen, front = selFirstFronts(unique_fits, pareto_fronts, k)
    if len(chosen) < k:
        chosen += selCrowdingRoundRobin(front, unique_fits, crowding_dist, k - len(chosen))

    for ind in chosen:
        ind.fitness.crowding_dist = crowding_dist[ind.fitness.wvalues]
        ind.fitness.rank = ranks[ind.fitness.wvalues]

    return chosen

def selTournamentFitnessDCD(individuals, k):
    """procedure UFTournSelection
    """
    def tourn(ind1, ind2):
        if ind2.fitness.rank < ind1.fitness.rank:
            return ind2
        elif ind1.fitness.rank < ind2.fitness.rank:
            return ind1

        if ind1.fitness.crowding_dist < ind2.fitness.crowding_dist:
            return ind2
        elif ind1.fitness.crowding_dist > ind2.fitness.crowding_dist:
            return ind1

        return choice((ind1, ind2))    

    unique_fits = defaultdict(list)
    for ind in individuals:
        unique_fits[ind.fitness.wvalues].append(ind)
    fitnesses = unique_fits.keys()

    if len(fitnesses) == 1:
        return sample(individuals, k)

    chosen = []

    while len(chosen) < k:
        samples = sample(fitnesses, min((k-len(chosen))*2, len(fitnesses)))
        for fit1, fit2 in zip(samples[::2], samples[1::2]):
            ind1 = choice(unique_fits[fit1])
            ind2 = choice(unique_fits[fit2])
            chosen.append(tourn(ind1, ind2))

    return chosen

def isDominated(wvalues1, wvalues2):
    """Returns whether or not *wvalues1* dominates *wvalues2*.

    :param wvalues1: The weighted fitness values that would be dominated.
    :param wvalues2: The weighted fitness values of the dominant.
    :returns: :obj:`True` if wvalues2 dominates wvalues1, :obj:`False`
              otherwise.
    """
    not_equal = False
    for self_wvalue, other_wvalue in zip(wvalues1, wvalues2):
        if self_wvalue > other_wvalue:
            return False
        elif self_wvalue < other_wvalue:
            not_equal = True
    return not_equal

def sortFastND(fitnesses, k, first_front_only=False):
    """Sort the first *k* *fitnesses* according the the fast non-dominated
    sorting algorithm.

    :param fitnesses: A list of unique fitnesses to sort in front.
    :param k: The number of individuals to select.
    :param first_front_only: If :obj:`True` sort only the first front and
                             exit.
    :returns: A list of Pareto fronts (lists), with the first list being the
              true Pareto front.
    """
    if k == 0:
        return []

    N = len(fitnesses)
    pareto_fronts = []

    pareto_fronts.append([])
    pareto_sorted = 0
    dominating_fits = [0] * N
    dominated_fits = [list() for i in range(N)]

    # Rank first Pareto front
    for i, fit_i in enumerate(fitnesses):
        for j, fit_j in enumerate(fitnesses[i+1:], i+1):
            if isDominated(fit_i, fit_j):
                dominating_fits[i] += 1
                dominated_fits[j].append(i)
            elif isDominated(fit_j, fit_i):
                dominating_fits[j] += 1
                dominated_fits[i].append(j)
        if dominating_fits[i] == 0:
            pareto_fronts[-1].append(i)
            pareto_sorted += 1

    # Rank the next front until all individuals are sorted or the given
    # number of individual are sorted
    if not first_front_only:
        k = min(N, k)
        while pareto_sorted < k:
            pareto_fronts.append([])
            for indice_p in pareto_fronts[-2]:
                for indice_d in dominated_fits[indice_p]:
                    dominating_fits[indice_d] -= 1
                    if dominating_fits[indice_d] == 0:
                        pareto_fronts[-1].append(indice_d)
                        pareto_sorted += 1

    ranks = {}
    for i, front in enumerate(pareto_fronts):
        for j, idx in enumerate(front):
            front[j] = fitnesses[idx]
            ranks[fitnesses[idx]] = i

    return pareto_fronts, ranks

def selFirstFronts(unique_fits, pareto_fronts, k):
    chosen = []
    for front in pareto_fronts:
        size = sum(len(unique_fits[fit]) for fit in front)
        if len(chosen) + size <= k:
            chosen.extend(chain(*(unique_fits[fit] for fit in front)))
        else:
            break
    return chosen, front

## Last front selection functions
def selCrowdingRoundRobin(front, unique_fits, crowding_dist, k):
    """procedure LastFrontSelection"""
    front.sort(key=crowding_dist.__getitem__, reverse=True)
    chosen = []
    for fitness in cycle(front):
        if len(unique_fits[fitness]) > 0:
            ind = unique_fits[fitness].pop()
            chosen.append(ind)
            if len(chosen) == k:
                break
    return chosen

def computeCrowdingDist(pareto_fronts):
    """Compute a crowding distance for each fitness on each front. The
    function retuns a dictonnary where the key is the
    fitness and the value the crowding distance.
    """
    nobj = len(pareto_fronts[0][0])
    distances = defaultdict(float)

    for front in pareto_fronts:
        for i in range(nobj):
            front.sort(key=itemgetter(i))
            distances[front[-1]] = float("inf")
            distances[front[0]] = float("inf")
            if front[-1][i] == front[0][i]:
                continue
            norm = float(front[-1][i] - front[0][i]) * nobj
            for prev, cur, nex in zip(front[:-2], front[1:-1], front[2:]):
                distances[cur] += (nex[i] - prev[i]) / norm
    return distances

__all__ = ['selNSGA2', 'selTournamentFitnessDCD']

