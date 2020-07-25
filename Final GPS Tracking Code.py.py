#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
Provides functions for obtaining current
satellite coordinates
"""
from matplotlib import pyplot as plt
import time
import string
import pynmea2
import serial
import math
from datetime import datetime
import requests

def get_current_location(connection_port):
    """
    Returns a dictionary with the current location of 
    Sat given the port to which the gps tracker is 
    connected to the onboard pc.
    """
    my_coordinates = {}

    while True:
        port = connection_port
        ser = serial.Serial(port, baudrate=9600, timeout=0.5)
        data_out = pynmea2.NMEAStreamReader()
        new_data = ser.readline()

        if new_data[0:6] == "$GPRMC":
            new_msg = pynmea2.parse(new_data)
            my_coordinates["lat"] = new_msg.latitude
            my_coordinates["long"] = new_msg.longitude
            return my_coordinates


"""
Provides functions and equations for obtaining the
weather details of the current location.
They include the temperature, humidty, pressure,
wind speed and direction, and the air density at 
a given/obtained temperature
"""
constants = {
    'specific_gas_constant_dry_air': 287.058,
    'specific_gas_constant_water_vapor': 461.495
}


def get_current_weather(latitude, longitude):
    """
    Gets Weather of current location

    Given the coordinates (latitude, longitude)
    returns a dictionary of the atmospheric data and location

    Temperature is returned in kelvin
    WIndspeed is in metres/second
    """
    try:
        weather_data = {}
        speed_dict = {}
        url = f'https://api.openweathermap.org/data/2.5/onecall?lat={latitude}&lon={longitude}&appid=c22f851ba43db08fb5f4938a532fa010'
        res = requests.get(url).json()
        weather_data['location'] = res['timezone']
        weather_data['date_time'] = str(datetime.fromtimestamp(res['current']['dt']))  # date and time, in local time
        weather_data['air_temperature'] = res['current']['temp']  # kelvin
        weather_data['air_pressure'] = res['current']['pressure']
        weather_data['humidity'] = res['current']['humidity']
        speed_dict['speed'] = res['current']['wind_speed']  # metres/second
        speed_dict['degree'] = res['current']['wind_deg']
        weather_data['wind_speed'] = speed_dict
        return weather_data

    except requests.ConnectionError:
        return print('Kindly check your intenet connection and try again')


# print(get_current_weather(0.32, 32.58))

def kelvin_to_celisius(kelvin_value):
    """
    Converts form Kelvins to degrees Celsius
    """
    return kelvin_value-273.15


def calculate_vapor_pressure(temperature, humidity):
    """
    Returns the vapor pressure given the
    temperature and humidity form weatehr API
    service
    """
    ratio = (7.5*kelvin_to_celisius(temperature))/(kelvin_to_celisius(temperature)+237.3)
    saturation_vapour = 6.1078*pow(10, ratio)
    return saturation_vapour*humidity


def calculate_dry_air_pressure(air_pressure, vapor_pressure):
    """
    Returns the dry air pressure given
    the total pressure from weather api
    service and the vapor pressure
    """
    return air_pressure - vapor_pressure

def calculate_air_density(vapor_pressure, dry_air_pressure, temperature):
    """
    Returns the air density which is a
    summation of the dry air and 
    air vapor densities
    """
    dry_air_part = (dry_air_pressure)/(constants['specific_gas_constant_dry_air']*kelvin_to_celisius(temperature))
    vapor_part = (vapor_pressure)/(constants['specific_gas_constant_water_vapor']*kelvin_to_celisius(temperature))
    return dry_air_part+vapor_part


"""
Provides functions used for obtaining the total mass
of the satellite given input data/parameters
"""
def calculate_additional_inertia_mass(air_density, balloon_volume):
    """
    Returns the additional mass developed as a result
    of inertia during motion of the satellite
    """

    return 0.5*air_density*balloon_volume


def calculate_helium_mass(helium_density, balloon_volume):
    """
    Returns the mass of helium per different volume of helium 
    """
    return helium_density*balloon_volume


def calculate_gross_mass(payload, film, helium, inertia):
    """
    Returns the total mass of the whole satellite body
    """
    return payload+film+helium+inertia


"""
Provides functions used in calculating the
aerodynamic forces expreinced during the flight
"""
def calculate_lift(balloon_volume, helium_density, air_density):
    """
    Returns the lift of the satellite at different
    stages/locations during the flight
    """
    # returns the negative of the equation because
    # air_density is more than helium_density
    # this helps us attain a positive value to be
    # used durong plotting purposes
    return abs(balloon_volume*(helium_density-air_density)*9.81)


def calculate_ascension_rate(balloon_volume, air_density, helium_density, drag_coefficient, surface_area):
    """
    Returns the vertical speed of balloon during flight
    """
    #included a 2 bse it had been missed and absolute ascension rate is calculated.
    return abs(math.sqrt(abs((2*balloon_volume*(air_density-helium_density)*9.81)/(drag_coefficient*air_density*surface_area))))


def calculate_relative_velocity(ascension_rate, wind_speed, wind_direction):
    """
    Returns the resultant velocity of balloon
    """
    v_x = wind_speed*math.cos(wind_direction)
    v_y = wind_speed*math.sin(wind_direction)
    return math.sqrt(pow(v_x, 2)+pow(v_y, 2)+pow(ascension_rate, 2))


def calculate_drag(drag_coefficient, air_density, relative_velocity, surface_area):
    """
    Returns the drag forces experienced by the balloon during flight
    """
    return 0.5*drag_coefficient*air_density*surface_area*pow(relative_velocity, 2)


def calculate_weight(gross_mass):
    """
    Returns the weight of the whole satellite body
    """
    return gross_mass*9.81


def calculate_new_volume(prev_pressure, prev_temp, prev_volume, new_pressure, new_temp):
    """
    Returns the current volume when there is 
    a change in pressure and temperature based on
    the ideal gas law
    """
    return (prev_pressure*prev_volume*new_temp)/(prev_temp*new_pressure)

#added surface area function since change in volume results in change in surface area.
def calculate_new_surface_area(new_volume):
    """
    Returns the current surface area when there is
    change in pressure and temperature
    """
    radius = ((3*new_volume)/(4*3.14))**(1/3)
    return 3.14*radius

"""
main code as a result of integration of the other 
sub codes
"""

helium_density = 0.1785
drag_coefficient = 0.285
balloon_film_mass = float(input('Input the film mass\n'))
payload_mass = float(input('Input the payload mass\n'))
# helium_mass = float(input('Input the helium mass\n'))
balloon_volume = float(input('Input the balloon volume\n'))
# balloon_diameter = float(input('Input the balloon diameter\n'))
balloon_surface_area = float(input('Input the balloon surface area\n'))
#drag_coefficient = float(input('Input the drag co-efficient\n'))

#for use in real world for coolecting coordinate data from the dictionary to a list
#uncomment with real world
#coordinates = list(my_coordinates.item())
#print(coordinates)
coordinates = [[0.165903, 32.438140], [0.366667, 32.440493], [0.621925, 32.310035], [0.901453, 32.127055],
              [1.212014, 31.811922], [1.314210, 31.599010]]

pressure_values = []
temperature_values = []
volume_values = []

# for plotting
lift_values = []
ascension_rate_values = []
relative_velocity_values = []
drag_values = []
weight_values = []

 # changed the name of the function used in the earlier code with calculate
for location in coordinates:
    if coordinates.index(location) == 0:
        weather = get_current_weather(location[0], location[1])
        pressure_values.append(weather['air_pressure'])
        temperature_values.append(weather['air_temperature'])
        volume_values.append(balloon_volume)
        vapor_pressure = calculate_vapor_pressure(weather['air_temperature'], weather['humidity'])
        dry_air_pressure = calculate_dry_air_pressure(weather['air_pressure'], vapor_pressure)
        air_density = calculate_air_density(vapor_pressure, dry_air_pressure, weather['air_temperature'])
        lift = calculate_lift(balloon_volume, helium_density, air_density)
        ascension_rate = calculate_ascension_rate(balloon_volume, air_density, helium_density, drag_coefficient, balloon_surface_area)
        relative_velocity = calculate_relative_velocity(ascension_rate, weather['wind_speed']['speed'], weather['wind_speed']['degree'])
        drag = calculate_drag(drag_coefficient, air_density, relative_velocity, balloon_surface_area)
        additional_inertia_mass = calculate_additional_inertia_mass(air_density, balloon_volume)
        helium_mass = calculate_helium_mass(helium_density, balloon_volume)
        gross_mass = calculate_gross_mass(payload_mass, balloon_film_mass, helium_mass, additional_inertia_mass)
        weight = calculate_weight(gross_mass)
        lift_values.append(round(lift, 4))
        ascension_rate_values.append(round(ascension_rate, 4))
        relative_velocity_values.append(round(relative_velocity, 4))
        drag_values.append(round(drag, 4))
        weight_values.append(round(weight, 4))
    else:
        # account for change in volume due to temperature and pressure changes
        weather = get_current_weather(location[0], location[1])
        pressure_values.append(weather['air_pressure'])
        temperature_values.append(weather['air_temperature'])
        vapor_pressure = calculate_vapor_pressure(weather['air_temperature'], weather['humidity'])
        dry_air_pressure = calculate_dry_air_pressure(weather['air_pressure'], vapor_pressure)
        air_density = calculate_air_density(vapor_pressure, dry_air_pressure, weather['air_temperature'])
        new_volume = calculate_new_volume(pressure_values[coordinates.index(location)-1], temperature_values[coordinates.index(location)-1], volume_values[coordinates.index(location)-1], weather['air_pressure'], weather['air_temperature'])
        volume_values.append(round(new_volume, 4))
        lift = calculate_lift(new_volume, helium_density, air_density)
        new_surface_area = calculate_new_surface_area(new_volume)
        ascension_rate = calculate_ascension_rate(new_volume, air_density, helium_density, drag_coefficient, new_surface_area)
        relative_velocity = calculate_relative_velocity(ascension_rate, weather['wind_speed']['speed'], weather['wind_speed']['degree'])
        drag = calculate_drag(drag_coefficient, air_density, relative_velocity, new_surface_area)
        inertia_mass = calculate_additional_inertia_mass(air_density, new_volume)
        helium_mass = calculate_helium_mass(helium_density, new_volume)
        gross_mass = calculate_gross_mass(payload_mass, balloon_film_mass, helium_mass, inertia_mass)
        weight = calculate_weight(gross_mass)
        lift_values.append(round(lift, 4))
        ascension_rate_values.append(round(ascension_rate, 4))
        relative_velocity_values.append(round(relative_velocity, 4))
        drag_values.append(round(drag, 4))
        weight_values.append(round(weight, 4))

print(f'Lift values: {lift_values}')
print(f'Ascension values: {ascension_rate_values}')
print(f'Velocity values: {relative_velocity_values}')
print(f'Drag values: {drag_values}')
print(f'Weight values: {weight_values}')
print(f'Pressure values: {pressure_values}')
print(f'Temperature values: {temperature_values}')
print(f'Volume values: {volume_values}')

"""
Graphs showing the results obtained

!Note: the actual graphs were obtained using jupyter notebooks due to vscode setup challenges
Uncomment the sections below and run
"""

locations = [1, 2, 3, 4, 5, 6]
# lift per given location
fig, ax = plt.subplots()
ax.plot(locations, lift_values)
ax.set_xlabel('Location')
ax.set_ylabel('Lift')

# ascension rate per given location
fig, ax = plt.subplots()
ax.plot(locations, ascension_rate_values)
ax.set_xlabel('Location')
ax.set_ylabel('Ascension')

# relative velocity per given location
fig, ax = plt.subplots()
ax.plot(locations, relative_velocity_values)
ax.set_xlabel('Location')
ax.set_ylabel('Relative velocity')

# volume changes with temperature changes
fig, ax = plt.subplots()
ax.plot(temperature_values, volume_values)
ax.set_xlabel('Temperature')
ax.set_ylabel('Volume')

# volume changes with pressure changes
fig, ax = plt.subplots()
ax.plot(pressure_values, volume_values)
ax.set_xlabel('Pressure')
ax.set_ylabel('Volume')

# ascension rate changes with temperature changes
fig, ax = plt.subplots()
ax.plot(temperature_values, ascension_rate_values)
ax.set_xlabel('Temperature')
ax.set_ylabel('Ascension rate')


# In[ ]:





# In[ ]:




