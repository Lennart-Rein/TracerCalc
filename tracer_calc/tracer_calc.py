# ---------------------------------------------------------
# Coding: UTF-8
# ---------------------------------------------------------
#
# License: Apache 2.0
# 
# Author: Lennart Rein
#
# This module provides methods for planning tracer tests.
# ---------------------------------------------------------

import math
import matplotlib.pyplot as plt

class TimeUnit:
    """
    Definitions of time units. Use the parameters to specify the output unit. 
    Use the methods to turn any time unit to seconds.
    """
    MILLISECONDS = "Milliseconds"
    SECONDS = "Seconds"
    MINUTES = "Minutes"
    HOURS = "Hours"
    DAYS = "Days"
    YEARS = "Years"

    def milliseconds(self, count):
        """
        Converts the input milliseconds to seconds
        """
        return count / 1000

    def seconds(self, count):
        """
        Converts the input seconds to seconds. This method exists for readability improvements.
        """
        return count

    def minutes(self, count):
        """
        Converts the input minutes to seconds
        """
        return count * 60

    def hours(self, count):
        """
        Converts the input hours to seconds
        """
        return count * 3600

    def days(self, count):
        """
        Converts the input days to seconds
        """
        return count * 86400

    def years(self, count):
        """
        Converts the input years to seconds
        """
        return count * 365 * 86400

class MassUnit:
    """
    Definitions of mass units. Use the parameters to specify the ouput unit.
    Use methods to turn any mass unit into milligrams.
    """
    MILLIGRAMS = 1
    GRAMS = 2
    KILOGRAMS = 3

    def milligrams(self, count):
        """
        Converts the input milligrams to milligrams. This method exists for readability improvements.
        """
        return count

    def grams(self, count):
        """
        Converts the input grams to milligrams.
        """
        return count*1000

    def kilograms(self, count):
        """
        Converts the input kilograms to milligrams.
        """
        return count * 1000000

class LengthUnit:
    """
    Definitions of distance units. Use the methods to turn any distance unit 
    into inches.
    """

    def cm(count):
        """
        Converts the input centimeters to inches.
        """
        return count/2.58

    def inch(count):
        """
        Converts the input inches to inches. This method exists for readability improvements.
        """
        return count
    

def plot_models(
    breakthough_model: tuple[list[float], list[tuple[float, list[float]]]],
    graph_time_unit: TimeUnit,
    plot_width: float = LengthUnit.cm(17),
    plot_height: float = LengthUnit.cm(10),
    font: str = "Arial",
    title: str = None,
    title_size: float = 11,
    legend_size: float = 9,
    axis_label_size: float = 11,
    tick_label_size: float = 9
):
    if len(breakthough_model[0]) == 0 or len(breakthough_model[0]) == 0:
        ValueError("Model is empty. Try calculating the model using 'model_breakthrough_curve_1D' or 'model_breakthrough_curve_2D' beforehand.")
    if graph_time_unit == TimeUnit.YEARS:
        x_labels = [time / (86400*365) for time in breakthough_model[0]]
    elif graph_time_unit == TimeUnit.DAYS:
        x_labels = [time / (86400) for time in breakthough_model[0]]
    elif graph_time_unit == TimeUnit.HOURS:
        x_labels = [time / 3600 for time in breakthough_model[0]]
    elif graph_time_unit == TimeUnit.MINUTES:
        x_labels = [time / 60 for time in breakthough_model[0]]
    elif graph_time_unit == TimeUnit.SECONDS:
        x_labels = breakthough_model[0]
    elif graph_time_unit == TimeUnit.MILLISECONDS:
        x_labels = [time * 1000 for time in breakthough_model[0]]
    else:
        raise ValueError("graph_time_unit must be any of the constants of TimeUnit")
    
    line_styles = []
    for i in range(len(breakthough_model[1])-1):
        line_styles.append((0, (2 * i + 1, 1)))
    line_styles.append("-")

    plt.figure(figsize=(plot_width, plot_height), dpi=300)
    plt.plot(x_labels, breakthough_model[1][0][1], marker=None, linestyle=None, color="1", label="Peclet Numbers")
    for index, value_tuple in enumerate(breakthough_model[1]):
        plt.plot(x_labels, value_tuple[1], marker=None, linestyle=line_styles[index], color="0", label=f"{value_tuple[0]}")

    plt.xlabel(graph_time_unit, fontname=font, fontsize=axis_label_size)
    plt.ylabel("Tracer concentration in mg/m$^3$", fontname=font, fontsize=axis_label_size)
    plt.xticks(fontname=font, fontsize=tick_label_size)
    plt.yticks(fontname=font, fontsize=tick_label_size)

    plt.title(title, fontname=font, fontsize=title_size)
    plt.grid(True)
    plt.legend(loc="upper right", prop={"family": font, "size": legend_size}, fancybox=False)
    plt.savefig(r"C:\Users\lenna\Documents\Uni\TracerTech\code\test.png")
    plt.show()


def model_breakthrough_curve_1D(
    start_time: float,
    end_time: float,
    timestep_count: int,
    distance: float,
    distance_velocity: float,
    flow_rate: float,
    tracer_mass: float,
    peclet: list[float],
) -> tuple[list[float], list[tuple[float, list[float]]]]:
    """
    This method models the breakthrough curve of a tracer test using the 1D 
    approximation of the Advection Dispersion Equation (ADE).

    :param start_time: 
        time the model should start in seconds.
    :param end_time: 
        time the model should end in seconds.
    :param timestep_count: 
        The amount of timesteps to be calculated inbetween the start and end time. 
        Higher numbers produce a smoother curve at the cost of performance.
    :param distance:
        distance for which the ADE should be solved. This distance will be used to calculate
        the dispersivity in meters.
    :param distance_velocity:
        distance velocity of the carrier medium of the tracer in m/s.
    :param discharge:
        discharge of the carrier medium of the tracer in m³/s.
    :param tracer_mass:
        mass of the tracer in milligrams.
    :param peclet:
        list of peclet numbers for each of which a break through curve should be modeled.

    :return: 
        a tuple with a list of the calculeted timesteps in seconds, as well as a list 
        of all the lists of calculated tracer concentrations. (timesteps, list[list[concentrations]]).
    """
    y_arrays = []
    delta = end_time - start_time
    timestep_size = delta / timestep_count
    x_values = [timestep * timestep_size for timestep in range(timestep_count)]
    for pecl in peclet:
        dispersivity = pecl * distance
        dispersion_coefficient = dispersivity * distance_velocity

        y_values = []

        for time in x_values:
            factor1 = tracer_mass / flow_rate
            factor2 = distance / math.sqrt(4 * math.pi * dispersion_coefficient * math.pow(max(time, 0.00001),3))
            factor3 = math.exp(-math.pow(distance - distance_velocity * max(time, 0.00001), 2)/(4 * dispersion_coefficient * max(time, 0.00001)))
            y_values.append(
                factor1 * factor2 * factor3
            )
        
        y_arrays.append(
            (pecl, y_values)
        )

    return (x_values, y_arrays)


def model_breakthrough_curve_2D(
    start_time: float,
    end_time: float,
    timestep_count: int,
    distance: float,
    aquifer_porosity: float,
    aquifer_thickness: float,
    distance_velocity: float,
    peclet: list[float],
    tracer_mass: float,
    width: float = None
) -> tuple[list[float], list[tuple[float, list[float]]]]:
    """
    This method models the breakthrough curve of a tracer test in an aquifer 
    using the 2D approximation of the Advection Dispersion Equation (ADE).

    :param start_time: 
        time the model should start in seconds.
    :param end_time: 
        time the model should end in seconds.
    :param timestep_count: 
        The amount of timesteps to be calculated inbetween the start and end time. 
        Higher numbers produce a smoother curve at the cost of performance.
    :param distance:
        distance for which the ADE should be solved. This distance will be used to calculate
        the dispersivity in meters.
    :param aquifer_porosity:
        the porosity of the aquifer in meters.
    :param aquifer_thickness:
        the thickness of the aquifer in meters.
    :param peclet:
        list of peclet numbers for each of which a break through curve should be modeled.
    :param tracer_mass:
        the mass of tracer used in milligrams.
    :param plume_width:
        the average width of the tracer plume in meters. If no value or None is passed a width of 0.1 * distance
        will be used.

    :return: 
        a tuple with a list of the calculeted timesteps in seconds, as well as a list 
        of all the lists of calculated tracer concentrations. (timesteps, list[list[concentrations]]).
    """
    if width is None:
        width = distance/10

    value_arrays = []
    time_delta = end_time - start_time
    time_step_size = time_delta / timestep_count
    time_values = [timestep * time_step_size for timestep in range(timestep_count)]

    for pecl in peclet:
        dispersion = distance * pecl
        dispersion_lateral = distance_velocity * dispersion
        dispersion_transversal = dispersion_lateral * 0.1

        y_values = []
        for time in time_values:
            time = max(time, 0.000001)
            factor_1 = tracer_mass / (aquifer_porosity * aquifer_thickness)
            factor_2 = distance / (4 * math.pi * distance_velocity * math.pow(time, 2) * math.sqrt(dispersion_lateral * dispersion_transversal))
            factor_3 = math.exp(-math.pow(distance - distance_velocity * time, 2)/(4 * dispersion_lateral * time)-math.pow(width, 2)/(4 * dispersion_transversal * time))
        
            y_values.append(
                factor_1 * factor_2 * factor_3
            )

        value_arrays.append(
            (pecl, y_values)
        )
    
    return (time_values, value_arrays)


def estimate_volume_simple(
    distance: float,
    aquifer_thickness: float,
    aquifer_porosity: float
):
    """
    Simple, physically based method of estimating the affected groundwater volume
    of a tracer test.

    :param distance:
        distance between injection and abstraction in meters.
    :param aquifer_thickness:
        thickness of the aquifer in meters.
    :param aquifer_porosity:
        porosity of the aquifer.
    """
    return distance * distance/10 * aquifer_thickness * 0.5 * aquifer_porosity


def estimate_volume_well_test(
    distance: float,
    aquifer_thickness: float,
    aquifer_porosity: float,
    is_dipole: bool
):
    """
    Empirical formula for estimating the affected groundwater volume
    of a tracer test at an injection (and abstraction) well.

    :param distance:
        distance between injection and abstraction in meters.
    :param aquifer_thickness:
        thickness of the aquifer in meters.
    :param aquifer_porosity:
        porosity of the aquifer.
    """
    p = 1
    if is_dipole:
        p = 3

    return p * math.pi * distance * distance * aquifer_thickness * aquifer_porosity


def calculate_distance_velocity(
      gradient: float,
      hydraulic_conductivity: float,
      porosity: float  
):
    """
    Method for calculating the distance velocity.

    :param gradient:
        gradient of the hydraulic head in m/m.
    :param hydraulic_conductivity:
        hydraulic conductivity of the aquifer in m/s.
    :param porosity:
        porosity of the aquifer.
    """
    return gradient * hydraulic_conductivity / porosity


def estimate_aquifer_flow_rate(
    distance_velocity: float,
    distance: float,
    thickness: float
):
    """
    Method for estimating the average flow rate through an aquifer.

    :param distance_velocity:
        distance velocity of the water in m/s.
    :param distance:
        distance between injection and measurement in meters.
    :param thickness:
        thickness of the aquifer in meters.
    """
    # The cross sectional area is halfed because the tracer 
    # plume is though of as a three dimensional triangle.
    cross_section = distance * 0.1 * thickness * 0.5
    return distance_velocity * cross_section


def estimate_tracer_mass(
    volume: float,
    relative_fluorecense_yield: float,
    out_unit: MassUnit
):
    """
    This method estimates the tracer amount needed for a tracer test.

    :param volume:
        water volume the tracer will be in in m³.
    :param relative_fluorecense_yield:
        relative fluorecense to uranin as a float between 0 and 1.
    :param out_unit:
        unit the output should be in. Can be any constant of MassUnit.
    """
    a = 1 / relative_fluorecense_yield
    mg =  10 * 0.01 * a * volume

    if out_unit == MassUnit.MILLIGRAMS:
        return mg
    elif out_unit == MassUnit.GRAMS:
        return mg/1000
    elif out_unit == MassUnit.KILOGRAMS:
        return mg/1000000
    else:
        ValueError("out_unit must be any of the constants of MassUnit.")


def test():
    print("test")

def sample_model(
    model, 
    sample_times: list[float]
):  
    concentration_arrays = []
    for peclet_and_concentrations in model[1]:
        peclet = peclet_and_concentrations[0]
        concentrations = peclet_and_concentrations[1]

        sampled_concentrations = []
        for t in sample_times:
            for i in range(len(model[0]) - 1):
                if model[0][i] <= t <= model[0][i + 1]:
                    t1, c1 = model[0][i], concentrations[i]
                    t2, c2 = model[0][i + 1], concentrations[i + 1]
                    
                    y = c1 + (c2 - c1) * (t - t1) / (t2-t1)
                    sampled_concentrations.append(y)
                    break
        
        concentration_arrays.append(
            (peclet, sampled_concentrations)
        )

    return (sample_times, concentration_arrays)
    