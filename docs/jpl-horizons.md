
I have implemented part of the `astroquery` module to easily make simulations of objects on the JPL Horizons System, using `GET`/`POST` URL encoded requests. I have restricted the inputs to standardise the output for use in this project, but you should be able to get any body listed on the JPL Horizons System fairly easily.

> [!NOTE]
> The content on this page is unlikely to change and is stable.

## `horizons_query()`
Returns attributes from a single search query, i.e, 'Sun'.
 ### Parameters
  *  `searchquery` (`str`) - object ID or identifiable name of object, ie, `'Sun'`.
  *  `observer` (`str`) - observer position ID. see JPL Horizons Manual for more info, default is `'0'` (Sun/Solar System Barycentric Centre). It is best to make this the reference object of your system.
 * `time` (`str`) - time to get data from. **Note: You should keep this the same for all objects in a system.** format is in MJY (YYYY-MM-DD).
 * `num_type` (`type`) - type of numerical to output information in, if needing high accuracy. Default is `float`.
 * `return_type` (`str`) how to output the result, either as a `dict` (`'dict'`), `Body` (`'body'`) or printing the result (`'print'`).

> [!TIP]
>In the case that the `searchquery` returns multiple unique objects, a list of objects should be outputted where you should enter the ID of your chosen object into the function instead.

 ### Returns
 * This function will return a `Body` instance with the attributes of the queried object by default.

## `horizons_batch()`
`horizons_batch` is a function that makes it easy to iterate over multiple objects using the same constants and return an iterable containing the result. Has the exact same arguments as `horizons_query()` except the input queries should be an iterable of `str` objects.

## Examples:

### Create an list of `Body` instances representing all the major bodies in the solar system.
```python
bodies = list(horizons_query(obj_id) for obj_id in (
                '10','199','299','399','499','599','699','799','899'))
```
> [!TIP]
>this could actually be made more concise:
>```python
>bodies = horizons_batch(['10','199','299','399','499','599','699','799','899'])
>```

>```python
>print(bodies)
>```
Returns:
>```
>Querying "10" @JPL Horizons System
>             ...
>Querying "899" @JPL Horizons System
>Getting data from JPL Horizons: 100%|█████████████████████████████████████████████████████████████████████████████| 9/9 [00:00<00:00, 98.66queries/s]
>[Body("Sun (10)", m=1.9885e+30, r=(...),v=(...)), a=(0, 0, 0)),
>             ...
>Body("Neptune (899)", m=1.02409e+26, r=(...),v=(...)), a=(0, 0, 0))]
>```
### Printing attributes of a body as `Decimal` and Querying using text:
> running this: would give us an error, 
>```python
>from decimal import Decimal
>horizons_query('Pluto', num_type=Decimal, return_type='print')
>```
>Returns:
>```
>Querying "Pluto" @JPL Horizons System
>Traceback (most recent call last):
>             ...
>ValueError: Ambiguous target name; provide unique id:
>  ID#      Name                               Designation  IAU/aliases/other
>  -------  ---------------------------------- -----------  -------------------
>        9  Pluto Barycenter
>      999  Pluto                              134340
>```
> The error is given as the horizons server picked up two possible objects that it could query, and cannot work out which one you wanted.
>instead, after the error we can choose the correct ID for the body we want to return, i.e, 999.
>```python
>horizons_query('999', num_type=Decimal, return_type='print')
>```
>Returns:
>```
>Querying "999" @JPL Horizons System
>Object: Pluto (999)
>***********
>Mass: 13070000000000000000000.000 kg
>Radius: 1188300.0 m
>***********
>position: (2436011869867.81494140625, -4641552098041.7890625, -250706318802.340087890625) (m)
>velocity: (24682.81863152435471420176327228546142578125, -21445.24496392483706586062908172607421875, -1503.5218150639748273533768951892852783203125)(ms^-1)
>***********
>Query Date: 2023-11-03
>```

> [!IMPORTANT]
> If you were to choose ''Pluto Barycenter'`, this would raise an exception as the function would not be able to find the mass or radius as a barycenter is not a physical object. 
