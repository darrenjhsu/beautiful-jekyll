---
layout: post
title: Minimal distance between two points of two groups
subtitle: 5 implementations with increasing performance
tags: programming-note
comments: true
---



Somewhere in my research I encountered this interesting problem: I have a large group of atoms A, and another, smaller group of atoms B. I want to find, for each atom of A, the closest distance to any B. To add a bit of context, group A is water (6325 oxygen atoms), and group B is protein atoms (304 atoms). To add complexity, I have 5000 frames from a trajectory, where water can move but protein is fixed.

Back of the envelope calculation shows that there will be 5000 * 6325 * 304 = \~1e10 distance calculations. That is quite a lot, so it's definitely worth looking into how to make this calculation faster. I'll present several methods with increasing efficiency.

Let's set up the variables here: `Ppos` (protein position) is a 304-by-3 matrix that does not change, and `Wpos` (water position) is a 6325-by-3 matrix for every frame. To access water position for each frame, assume the following.

```python
for frame in range(5000):
  Wpos = getWaterPosition()

def getWaterPosition(size=6325):
  return np.random.rand(size,3)

def getProteinPosition(size=304):
  return np.random.rand(size,3)
```

Our goal is to accumulate the shortest distance into a matrix `distList`, and perform a histogram method to get the distribution.

### Approach 1: The really naive way

The really intuitive way would be to calculate, in loops, for each water position, how far they are from each protein position, and then see if it's closer than the previous value. The code read like this:

```python
import numpy as np
Ppos = getProteinPosition()
distList = []
for frame in range(5000):
  Wpos = getWaterPosition()
  for i in Wpos:
    minDistance = 5000 # It's quite long
    for j in Ppos:
      distance = np.sqrt(((Ppos[j] - Wpos[i])**2).sum(1))
      if distance < minDistance:  
        minDistance = distance
    distList.append(minDistance)
```

The speed of this code is about 14 s / frame. That's quite slow! Let's improve it a bit.


### Approach 2: Using np.linalg and broadcast

The next way is to calculate the distances of a water to all proteins. This can be done by using the broadcasting feature of `numpy`. 

```python
import numpy as np
import numpy.linalg as la
Ppos = getProteinPosition()
distList = []
for frame in range(5000):
  Wpos = getWaterPosition()
  for i in Wpos:
    distList.append(np.min(la.norm(i - Ppos, axis = 1)))
```

This code dramatically increased efficiency to 165 ms / frame. Not bad! But for a 5000 frame trajectory this will take about 14 minutes. Can we make it faster?

### Approach 3: Using more broadcasting

I found this method on [stack overflow](https://stackoverflow.com/questions/28687321).

Instead of calculating distances of one water to all protein atoms independently, we can use more broadcasting by creating a third dimension and transposing one of the matrices to have numpy calculate the distances of all water atoms to all protein atoms in once.  In this code, `np.sqrt(((Wpos[:, :, None] - Ppos[:, :, None].T) ** 2).sum(1)` gives a `len(Wpos)`-by-`len(Ppos)` matrix which is the distance for all pairs. Taking the min of that would give the minimal distance for each water atom to any protein atoms.

```python
import numpy as np
Ppos = getProteinPosition()
distList = []
for frame in range(5000):
  Wpos = getWaterPosition()
  distList.append(np.min(np.sqrt(
      ((Wpos[:, :, None] - Ppos[:, :, None].T) ** 2).sum(1)),axis=1)
    )
```

This code improves to 82 ms / frame, a 2x increase in speed. Good! But as you can see the code is becoming hard for people who're not fluent in linear algebra to understand. At this point let's explore other possibilities.

### Approach 4: Using numba to leverage C++ speed

I noticed that above functions only utilize one core on my computer, but it has 8! Can we find a way to use up all cores when doing the `numpy` calculation? Before getting to that, I found [this article](https://louisabraham.github.io/articles/broadcasting-and-numba.html) that shows the usage of `numba`. We revert to the method in Approach 2 but use the `@njit` decorator to have it compiled.

```python
import numpy as np
from numba import njit, prange
Ppos = getProteinPosition()
distList = []
for frame in range(5000):
  Wpos = getWaterPosition()
  distList.append(WPFunction(Wpos, Ppos))

@njit
def WPFunction(Wpos, Ppos):
  WPLen = len(Wpos)
  distList = np.zeros(WPLen)
  for i in prange(WPLen):
    dist = np.sqrt(((WPos - Ppos)**2).sum(1))
    distList[i] = np.min(dist)
  return distList
```

This method give a speed of 36 ms / frame (excluding compile time). It's already about 400x faster than Approach 1.

### Approach 5: Using numba to leverage C++ speed and parallelism

The above code still only uses one core. It's time to have it use all cores by adding the `(parallel=True)` argument.

```python
import numpy as np
from numba import njit, prange
Ppos = getProteinPosition()
distList = []
for frame in range(5000):
  Wpos = getWaterPosition()
  distList.append(WPFunction_parallel(Wpos, Ppos))

@njit(parallel=True)
def WPFunction_parallel(Wpos, Ppos):
  WPLen = len(Wpos)
  distList = np.zeros(WPLen)
  for i in prange(WPLen):
    dist = np.sqrt(((WPos - Ppos)**2).sum(1))
    distList[i] = np.min(dist)
  return distList
```

This code now runs at 6.6 ms / frame. Impressive! Is it possible that we make it even faster?

### Approach 5+1: Stack frames and process them all at once

As a final try, I realized that it is not necessary to process all frames separately as I'm only interested in the distribution. Therefore, I could stack all water positions together and have them processed as a huge matrix.

```python
import numpy as np
from numba import njit, prange
Ppos = getProteinPosition()
distList = []
stackFrames = 5000
stackCounter = 0
Wpos = getWaterPosition()
WLen = len(Wpos)
WposStack = np.zeros((WLen*stackFrames,3))
for frame in range(5000):
  Wpos = getWaterPosition()
  WposStack[stackCounter*WLen:(stackCounter+1)*WLen] = Wpos
  if (frame+1) % stackFrames == 0:
    distList.append(WPFunction(Wpos, Ppos))
    WposStack = np.zeros((WLen*stackFrames,3))

@njit(parallel=True)
def WPFunction_parallel(Wpos, Ppos):
  WPLen = len(Wpos)
  distList = np.zeros(WPLen)
  for i in prange(WPLen):
    dist = np.sqrt(((Wpos - Ppos)**2).sum(1))
    distList[i] = np.min(dist)
  return distList
```

This code runs at 5.9 ms / frame using all cores, excluding compile time. We've come all the way from 14 s / frame down to this, an increase of 2400x in speed! This shows how it is possible to optimize the code and achieve amazing efficiency.

Note: If the system size (number of water and protein atoms) grow, the difference in speed between Approach 5 and 5+1 becomes smaller as the former is starting to use all CPU resources for that increasingly large calculation.

Note 2: Another way of doing this problem is to use a kd-tree (`scipy.spatial.cKDTree`) whose speed is similar to Approach 5. However, my original goal is to find the closest distance to van der Waals (vdW) surfaces (some element-specific radii from the atom center), so kd-tree is not a good way to do so. In these approaches, we can calculate that distance by subtracting vdW radii from the `dist`.

[back](./)
