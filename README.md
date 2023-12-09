
<h1 style="text-align: center;"><code>nBody</code> Scientific Modelling Project alpha-v2.0</h1>

  By Will Roscoe -- 9th Dec 2023

> TODO
> - [ ] Tests for body and Engine including    
>   - [ ] internal testing for collisions
>   - [ ] internal testing for gravity
>   - [ ] body update
> - [ ] Update functions in line with new mpmath usages in
>    - [ ] body.py
>    - [ ] engine.py
>    - [ ] visual.py
> - [ ] Importing csvs per body (opt)
----

The goal of this project is to efficiently simulate the gravitational and kinetic interactions between physical objects and display them in a user friendly interface. The usage of this project is Object Oriented and mean't to be flexible and easy to use.

>### Required Packages:
>
>`math`, `re`, `datetime`, `decimal`, `multiprocessing` (Built in)
>
>`tqdm`, `numpy`, `matplotlib`, `cycler`, `astroquery`, `astropy`
>
> ---
> `scipy` is recommended but optional.


# Usage Guide 
There are 3 main objects which a user would interact with while creating a simulation, these are:
*  ### [**`Body`**](#body) : Representation of a physical object
*  ### [**`PhysEngine`**](#physeng) : Representation of a system of `Body` objects and comuputes the physics.
*  ### [**`Simulation`**](#simul) : Controls the output and the runtime of a `PhysEngine` object.
----
### About Objects and Types
* There are also several other classes, `Variable`, `HistoricVariable`, `Vector` and `Hisoric Vector` which are classed as having types `VarType` and `VectorType` respectively.

* `NoneType` is defined as the type of `None`

<h2 style="text-align: left;", custom-id='body'>Creating <code>Body</code> Instances</h2>

The `Body` class is the representation of a physical body, like a planet or star or particle. To create a `Body`, you must pass it some parameters:
```python
new_body = Body(mass: Numeric,
                init_pos: list | tuple,
                init_vel: list | tuple=(0,0,0),
                radius: Numeric = 0,
                bounce: Numeric = 0.999,
                color: str= None,
                identity:str = None):
```

 - `init_pos` and `init_vel` must be passed a **list or tuple or Vector Object** with 3 values representing the cartesian components x,y,z in **m** and **ms^-1**.

 - `mass` is required and should be passed as a **numeric**, representing the total mass of the body in **kg**.

> - `radius` is an optional parameter, and if left as `None`, the object will not have collision dectection with other bodies or physical bounds. The units for `radius` is **m**.
> - `identity` is also optional, and defines what the label is for the body in the output. It must be passed a **string**. If left as `None`, a placeholder name will be created.
>


<h2 style="text-align: left;",custom-id='physeng'>The <code>PhysEngine</code> Instance</h2>

A `PhysEngine` instance computes and evaluates the attibutes of a set of bodies, and computes the next attributes taking into account the rest of the bodies in the instance.

To create a `PhysEngine` you can pass it a time interval, `dt`, or leave it as default.
```python
phys = PhysEngine(dt: int | float = 1, 
                  checking_range = int = 3):
```
You must then load the bodies using the function: 
```python
phys.attach_bodies(new_bodies:list | tuple) 
    #you must pass the bodies in a list or tuple.
```
>Optionally, you can create some infinite planes parallel to a plane by passing:
>```python
>phys.create_plane(self, const_axis='z', const_val = 0)
>```
>or induce a global accelerations by passing:
>```python
>phys.create_acceleration(accel_vector: tuple | list | VectorType)
>```
>or make bodies relative to a body by passing:
>```python
>phys.make_relative_to(target_body: Body)
>```
At this point, the `PhysEngine` instance is fully prepared to compute the trajectories of the bodies it was passed, which can be done using `phys.evaluate()` for as many steps as neccesary, however it works best to use the `Simulation` class to output the product to a graphical interface.

<h2 style="text-align: left;", custom-id='simul'>Creating and Using a <code>Simulation</code> Instance</h2>

A `Simulation` instance enables a user to run a `PhysEngine` easier and output an interactive GUI of the bodies' trajectories through `matplotlib`. There is only one required argument to pass to a `Simulation` instance, which is the `PhysEngine` containing the `Body` instances.
```python
sim = Simulation(name: str = 'Nbody Simulation',
                engine: PhysEngine|NoneType = None,#REQUIRED
                focus_body: Body|NoneType = None,
                focus_range: int|float = 0.5,
                autoscale: bool = True,
                show_grid: bool = True, 
                show_shadows: bool = False,
                show_acceleration: bool = False,
                show_velocity: bool = False,
                vector_size: int|float = 1,
                labelling_type: str = 'legend',
                body_model: str = 'dots',
                guistyle: str = 'default',
                ) -> None:
```
> The optional arguments are as follows:
>   | Argument    | Type  | Default        | Description |
>   | ----------- | ---- |----------- | ----------- |
>   | `name`| `str` | `'Nbody Simulation'`|Name of the Simulation or set of Bodies.|
>   | `focus_body`| `Body` or `NoneType` | `None`|Body to keep in centre of view.|
>   | `focus_range`| `int` or `float` | `0.5`| Distance to render to from location of `focus_body`.|
>   | `autoscale`| `bool` |`True`|whether or not to just autoscale axes, overrides most view options.|
>   | `show_grid`| `bool` |` True`|toggles grid visibility.|
>   | `show_shadows`|`bool` | `False`|toggles whether to plot a 2d projection on the xy plane.|
>   | `show_acceleration`| `bool` | `False`|toggles whether to draw acceleration vectors.|
>   | `show_velocity`| `bool` | `False`|toggles whether to draw velocity vectors.|
>   | `vector_size`| `int` or `float` | `1`|scalar multiplier of above vectors.|
>   | `labelling_type`| `str` | `'legend'`|either `'legend'` or `'label'`or `None`, defines whether to label the points, create a legend, or none.|
>   | `body_model`| `str` | `'dots'`|either `'dots'`, `'wireframe'` or `'surface'`, defines how to draw bodies with nonzero radii as a dot, or spherical surface.|
>   | `guistyle`| `str` | `'default'`|either `'default'` or `'dark'`, controls gui theming.|

After creating an instance, we can run the simulation and output the result using:
```python
sim.start(self,
            frames: int = None,
            interval: int|float = None,
            duration: int|float =None,
            fps: int|float=None):
```
you must pass two or more of the functions arguments, where `interval` is in seconds.

# The `horizons_query()` Function

I have implemented part of the `astroquery` module to easily make simulations of objects on the JPL Horizons System, using `GET`/`POST` URL encoded requests. I have restricted the inputs to standardise the output for use in this project, but you should be able to get any body listed on the JPL Horizons System fairly easily.

> ### Parameters
>  *  `searchquery` (`str`) - object ID or identifiable name of object, ie, `'Sun'`.
>  *  `observer` (`str`) - observer position ID. see JPL Horizons Manual for more info, default is `'0'` (Sun/Solar System Barycentric Centre). It is best to make this the reference object of your system.
> * `time` (`str`) - time to get data from. **Note: You should keep this the same for all objects in a system.** format is in MJY (YYYY-MM-DD).
> * `num_type` (`type`) - type of numerical to output information in, if needing high accuracy. Default is `float`.
> * `return_type` (`str`) how to output the result, either as a `dict` (`'dict'`), `Body` (`'body'`) or printing the result (`'print'`).

In the case that the `searchquery` returns multiple unique objects, a list of objects should be outputted where you should enter the ID of your chosen object into the function instead.

 ### Returns
 * This function will return a `Body` instance with the attributes of the queried object by default.

> ### `horizons_batch()`
> `horizons_batch` is a function that makes it easy to iterate over multiple objects using the same constants and return an iterable containing the result. Has the exact same arguments as `horizons_query()` except the input queries should be an iterable of `str` objects.

## Premade Examples
### `SolarSystemMB`: a Simulation of the Solar System's Major Bodies
You can initialise this simulation using its defaults by running the code below:

```python
solarsystem = SolarSystemMB()
solarsystem.start()
```
The options are the same as the `Simulation` Class, however autoscaling and focusing is locked to keep the Sun as the central body.

### `BouncingBalls`: Bouncing Balls in a Box, Colliding with each other and the walls.
Similarly, we initialise and run the simulation as below:
```python
bouncingballs = BouncingBalls()
bouncingballs.start()
```
This simulation is great at showcasing the ability for `Body` instances with non zero radius to collide with each other or physical planes. Note that the collisions may not always be accurate as the bodies collisions are modelled to be head on, and estimate collision points by evaluating positional information before and after the actual time of collision.
