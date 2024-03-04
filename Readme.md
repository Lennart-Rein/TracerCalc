# TracerCalc
  
1. About
2. Usage
3. Future Plans
4. Literature
5. License

### 1. About
TracerCalc is a light weight python module that simplifies the process of planning tracer tests. It provides simple methods to estimate aquifer parameters and needed tracer amounts, methods to model breakthrough curves according to 1D and 2D simplifications of the Advection-Dispersion-Equation (ADE), as well as methods to test different sampling strategies.  
The Module was developed as part of an exercise of the course "Tracer Techniques" at TU Darmstadt by Prof. Dr. Külls.

### 2. Usage

Aquifer Parameters can be estimated with easy to use utility methods.
```python
distance_velocity = calculate_distance_velocity(
    gradient=gradient, 
    hydraulic_conductivity=hydraulic_conductivity,
    porosity=aquifer_porosity
)

flow_rate = estimate_aquifer_flow_rate(
    distance_velocity=distance_velocity, 
    distance=distance, 
    thickness=aquifer_thickness
)

water_volume = estimate_volume_well_test(
    distance=distance,
    aquifer_thickness=aquifer_thickness,
    aquifer_porosity=aquifer_porosity,
    is_dipole=True
)

...
```

Needed tracer amounts can be estimated. Different tracers may be used by providing their fluorecense relative to that of uranine

```python
mass = estimate_tracer_mass(
    volume=water_volume,
    relative_fluorecense_yield=relative_fluorecense_yield,
    out_unit=MassUnit.GRAMS
)
```
Breakthrough curves can be modeled in using 1D and 2D approximations of the ADE.

```python
observation_distance = 50 

model_1d = model_breakthrough_curve_1D(
    start_time=0,
    end_time=TimeUnit().years(5),
    timestep_count=1000,
    distance=observation_distance,
    distance_velocity=distance_velocity,
    flow_rate=flow_rate,
    tracer_mass=MassUnit().grams(mass),
    peclet=[0.05, 0.1, 0.5, 1, 1.5, 2]
)

model_2d = model_breakthrough_curve_2D(
    start_time=0,
    end_time=TimeUnit().years(5),
    timestep_count=1000,
    distance=observation_distance,
    aquifer_porosity=aquifer_porosity,
    aquifer_thickness=aquifer_thickness,
    distance_velocity=distance_velocity,
    peclet = [0.05, 0.1, 0.5, 1, 1.5, 2],
    tracer_mass=MassUnit().grams(mass)
)
```

Furthermore the created models can be sampled in order to find valid sampling frequencies.

```python
sample_times = [(i+1) * TimeUnit().hours(12) for i in range(10)] # Sample 10 times every 12 hours
sample_times += [sample_times[-1] + (i+1) * TimeUnit().days(2) for i in range(10)] # then 10 times every 2 days

sample = sample_model(model=model_2d, sample_times=sample_times)
```

Models and sampled models can easily be plotted into beautiful graphs.

```python
plot_models(
    breakthough_model=model,
    graph_time_unit=TimeUnit.YEARS,
    title = "50 Meters, Sampled",
)
```

A jupyter notebook containing the codesnippets above and more can be found in the root folder of the repository.

### 3. Future Plans
- Better documentation of the used formulas
- Provide other commonly used formulas for mass estimation
- Better support for tracer tests in rivers

### 4. Literature
- Leibundgut, C., Maloszewski, P., Külls, C. (2009). Tracers in Hydrology. John Wiley & Sons, Ltd. DOI: 10.1002/9780470747148.

### 5. License
```
Copyright [2024] [Lennart Rein]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```


