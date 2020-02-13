# utility module for working with Apollo HFE data from various sources


import datetime as dt
import os
import copy
import re
import astropy.time as atime
import numpy as np
import pandas as pd


def load_hfe(pathname="."):
    data = {}

    # uncorrected 2005 NSSDC data
    for mission in ["a15", "a17"]:
        if not mission in data.keys():
            data[mission] = {}
        for probe in ["p1", "p2"]:
            data[mission][probe] = {}
            for sensor in [1, 2, 3, 4, 5]:
                filename = "{pathname}/{m}/{m}_hfe_{p}_{s}.tab".format(
                    pathname=pathname, m=mission, p=probe, s=sensor
                )
                columns = (
                    ["Time", "HTR", "TREF", "TC1", "TC2", "TC3", "TC4"]
                    if sensor == 3
                    else ["Time", "dT", "T"]
                )
                data[mission][probe][sensor] = pd.read_csv(
                    filename,
                    skiprows=2,
                    names=columns,
                    delim_whitespace=True,
                    skipinitialspace=True,
                )
    return data


def flag_missing_hfe(data):
    """Flag rows with any missing data values of -9999. Nagihara data don't use
    this convention (they simply don't include any missing points), but we
    create the flags column while we're at it.
    we do *not* exclude values where HTR data is missing; HTR data
    is essentially unusable, but thermocouple/reference bridge data isn't.
    """
    for mission in data.keys():
        for probe in data[mission].keys():
            for sensor in data[mission][probe].keys():
                flag_init = np.array(data[mission][probe][sensor]["Time"]).shape[0]
                data[mission][probe][sensor]["flags"] = pd.Series(
                    np.zeros(flag_init, dtype=np.int16),
                    index=data[mission][probe][sensor].index,
                )


                if sensor == 3:
                    data[mission][probe][sensor].loc[
                        (data[mission][probe][sensor]["Time"] == -9999)
                        | (data[mission][probe][sensor]["TC1"] == -9999)
                        | (data[mission][probe][sensor]["TC2"] == -9999)
                        | (data[mission][probe][sensor]["TC3"] == -9999)
                        | (data[mission][probe][sensor]["TC4"] == -9999),
                        "flags",
                    ] += 0b1
                else:
                    data[mission][probe][sensor].loc[
                        (data[mission][probe][sensor]["dT"] == -9999)
                        | (data[mission][probe][sensor]["T"] == -9999)
                        | (data[mission][probe][sensor]["Time"] == -9999),
                        "flags",
                    ] += 0b1


# flags time-inverted pairs of points in official data, then sorts all sets by time, which
# also fixes some misplaced time ranges in the Nagihara paper data. checks all sets,
# although this is a bit wasteful.


def manage_disordered_hfe(data):
    for mission in data.keys():
        for probe in data[mission].keys():
            for sensor in data[mission][probe].keys():
                if mission in ["a15", "a17"]:
                    flags = np.array(data[mission][probe][sensor]["flags"])
                    for i in np.arange(data[mission][probe][sensor]["Time"].size - 1):
                        if (
                            data[mission][probe][sensor]["Time"][i + 1]
                            - data[mission][probe][sensor]["Time"][i]
                            <= 0
                            and data[mission][probe][sensor]["Time"][i] != -9999
                            and data[mission][probe][sensor]["Time"][i + 1] != -9999
                        ):
                            flags[i] += 0b10000000000
                            flags[i + 1] += 0b10000000000
                    data[mission][probe][sensor]["flags"] = flags
                data[mission][probe][sensor] = (
                    data[mission][probe][sensor]
                    .sort_values(by="Time")
                    .reset_index(drop=True)
                )


# functions for writing out reduced sets.


def add_crlf_ending(filename):
    """Converts final terminator of file to CRLF.
    """
    if os.linesep != "\r\n":
        with open(filename, "rb") as file:
            needs_terminator = file.read()
        has_terminator = needs_terminator[:-1] + b"\r\n"
        with open(filename, "wb") as file:
            file.write(has_terminator)


def write_clean_hfe(data, outpath=".", version=""):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        os.mkdir(outpath + "/a15")
        os.mkdir(outpath + "/a17")
    data_clean_out = copy.deepcopy(data)
    for mission in data_clean_out.keys():
        for probe in data_clean_out[mission].keys():
            for sensor in data_clean_out[mission][probe].keys():
                if "dT_corr" in data_clean_out[mission][probe][sensor].columns:
                    data_clean_out[mission][probe][sensor] = data_clean_out[mission][
                        probe
                    ][sensor].reindex(columns=["Time", "T", "dT", "dT_corr", "flags"])
                elif sensor == 3:
                    data_clean_out[mission][probe][sensor] = data_clean_out[mission][
                        probe
                    ][sensor].reindex(
                        columns=[
                            "Time",
                            "HTR",
                            "TREF",
                            "TC1",
                            "TC2",
                            "TC3",
                            "TC4",
                            "flags",
                        ]
                    )
                else:
                    data_clean_out[mission][probe][sensor] = data_clean_out[mission][
                        probe
                    ][sensor].reindex(columns=["Time", "T", "dT", "flags"])

                # convert to the IBM 1130-equivalent number format, but retain millisecond precision
                # for Nagihara data

                if mission[3:4] == "_":
                    for column in data_clean_out[mission][probe][sensor].columns:
                        if column == "Time":
                            data_clean_out[mission][probe][sensor][column] = (
                                data_clean_out[mission][probe][sensor][column]
                                .apply("{:.11E}".format)
                                .str.pad(17, "right")
                            )
                        else:
                            data_clean_out[mission][probe][sensor][column] = (
                                data_clean_out[mission][probe][sensor][column]
                                .apply("{:.7E}".format)
                                .str.pad(14, "right")
                            )
                else:
                    for column in data_clean_out[mission][probe][sensor].columns:
                        data_clean_out[mission][probe][sensor][column] = (
                            data_clean_out[mission][probe][sensor][column]
                            .apply("{:.7E}".format)
                            .str.pad(14, "right")
                        )

                    # the following series of ugly regex substitutions is required because
                    # of deficiencies in fixed-width table formatting in pandas, exacerbated
                    # by a regression in pandas.to_string 0.25 that inserts leading
                    # spaces in non-indexed output.

                table = re.sub(
                    r"\n\s(\d|-)",
                    r"\n\1",
                    data_clean_out[mission][probe][sensor].to_string(
                        index=False, justify="left"
                    ),
                )
                table = re.sub(r" (?= (\d|-))", r"", table)
                # create correct line endings when script run from any major OS
                table = re.sub(r"(\n)|(\r\n)|(\r)", r"\r\n", table)
                
                filename = "{outpath}/{m}{p}f{s}{v}.tab".format(
                            outpath=outpath + "/" + mission[0:3],
                            m=mission,
                            p=probe,
                            s=sensor,
                            v=version,
                            )
                with open(
                    filename,
                    "w",
                ) as output_file:
                    print(table, file=output_file)
                add_crlf_ending(filename)

def write_split_hfe(data, outpath=".", version=""):
    data_split_out = copy.deepcopy(data)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for mission in data_split_out.keys():
        for probe in data_split_out[mission].keys():
            for sensor in data_split_out[mission][probe].keys():
                data_split_out[mission][probe][sensor]["Time"] = (
                    data_split_out[mission][probe][sensor]["Time"]
                    .dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
                    .str.slice(0, -3)
                    + "Z"
                )
                filename = "{outpath}/{m}{p}f{s}{v}_split.tab".format(
                    outpath=outpath, m=mission, p=probe, s=sensor, v=version
                )
                data_split_out[mission][probe][sensor].to_csv(
                    filename, index=False, line_terminator="\r\n",
                )

def write_deep_hfe(data, outpath=".", version=""):
    data_deep_out = copy.deepcopy(data)
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for mission in data_deep_out.keys():
        for probe in data_deep_out[mission].keys():
            data_deep_out[mission][probe]["Time"] = (
                data_deep_out[mission][probe]["Time"]
                .dt.strftime("%Y-%m-%dT%H:%M:%S.%f")
                .str.slice(0, -3)
                + "Z"
            )
            filename = "{outpath}/{m}{p}{v}_depth.tab".format(
                outpath=outpath, m=mission, p=probe, v=version
            )
            data_deep_out[mission][probe].to_csv(
                filename, index=False, line_terminator="\r\n",
            )

# Functions for interpreting data released by Nagihara et. al along with their 2018 paper
# "Examination of the Long-Term Subsurface Warming Observed at the Apollo 15 and 17 Sites
# Utilizing the Newly Restored Heat Flow Experiment Data From 1975 to 1977." Output intended
# primarily as intermediate data for further correction and reduction by other utilities in
# this module.

# The original reduced Apollo HFE data uses an epoch time format: milliseconds from
# 24 hours before the beginning of the mission's start year. This is December 31, 1970
# for Apollo 15, and December 31, 1971 for Apollo 17.
# These functions convert between Gregorian time and mission epoch time.

# epoch-to-Gregorian functions. placing in TAI to avoid leap second weirdness.
# excessive digits of precision are for parity with internal astropy values.


def a15_to_greg(time):
    epoch = dt.datetime(1970, 12, 31, 0, 0, 8, 943570)
    return epoch + dt.timedelta(milliseconds=time)


def a17_to_greg(time):
    epoch = dt.datetime(1971, 12, 31, 0, 0, 9, 889650)
    return epoch + dt.timedelta(milliseconds=time)


# Gregorian-to-epoch functions (assume TAI input)


def a15_time(time):
    if not time is None:
        epoch = dt.datetime(1970, 12, 31, 0, 0, 8, 943570)
        return (time - epoch).total_seconds() * 1000
    return None


def a17_time(time):
    if not time is None:
        epoch = dt.datetime(1971, 12, 31, 0, 0, 9, 889650)
        return (time - epoch).total_seconds() * 1000
    return None


# utility functions for converting between TAI and UTC. this allows us to use astropy
# to deal with leap seconds without later doing large table comparisons between astropy
# Time objects (much slower than using datetime).

# 9999 flags empty rows in nagihara 2018 spreadsheet.


def tai_to_utc(time):
    return atime.Time(time, scale="tai").utc.datetime


def utc_to_tai(time):
    if not time == dt.datetime(9999, 1, 1, 0, 0, 0):
        return atime.Time(time, scale="utc").tai.datetime


# silly utility functions for vectorizing over Nagihara datelists. 9999 flags empty rows.


def seconds_interval(time):
    if not np.isnan(time):
        return dt.timedelta(seconds=time)
    return dt.timedelta(0)


def days_since_year(day, year):
    if not np.isnan(day):
        return dt.datetime(year, 1, 1) + dt.timedelta(days=(int(day) - 1))
    return dt.datetime(9999, 1, 1, 0, 0, 0)


# Nagihara PDS release uses DOY format; this simply breaks it up to datetime equivalent.
# might be better to use strptime.


def nagihara_doy_to_dt(nagihara_time):
    year = dt.datetime(int(nagihara_time[0:4]), 1, 1, 0, 0)
    day_of_year = dt.timedelta(days=int(nagihara_time[5:8]) - 1)
    hour = dt.timedelta(hours=int(nagihara_time[9:11]))
    minute = dt.timedelta(minutes=int(nagihara_time[12:14]))
    second = dt.timedelta(seconds=float(nagihara_time[15:21]))
    return year + day_of_year + hour + minute + second


# main function for ingesting data from Nagihara et al. 2018


def ingest_nagihara_2018(
    nagihara_data=None, spreadsheet_path="./source/nagihara/jgre.xlsx"
):
    nagihara_datafiles = ["a15_1975", "a17_1975", "a17_1976", "a17_1977"]
    nagihara_spreadsheet = {}
    for file in enumerate(nagihara_datafiles):
        nagihara_spreadsheet[file[1]] = pd.read_excel(
            spreadsheet_path, file[0], header=0, delim_whitespace=True
        )

    # The 1975 data from this paper has been superseded by the related July 2019 PDS release.
    # We read in the whole spreadsheet for convenience, but only actually ingest
    # the 1976 and 1977 data.
    # as elsewhere, m is mission (here mission by year), p is probe,
    # s is 'file,' i.e. bridge identifier, which conveniently here is only ever '1.'

    nagihara_data = {}
    for mission in nagihara_spreadsheet.keys():
        if not mission[4:] == "1975":
            if not mission in nagihara_data.keys():
                nagihara_data[mission] = {}
            for probe in ["p1", "p2"]:
                if not probe in nagihara_data[mission].keys():
                    nagihara_data[mission][probe] = {}
                for sensor in [1]:
                    if not sensor in nagihara_data[mission][probe].keys():
                        nagihara_data[mission][probe][sensor] = {}

                    # The non-deprecated data from this paper is from the upper gradient bridges
                    # of both Apollo 17 probes.

                    # we initially read each 'file' (bridge) as a dict of arrays for convenience.

                    # Time

                    # the original HFE dataset reduced time data into an epoch time
                    # format as described above. Nagihara et al. reduced time into days
                    # from the beginning of the then-current calendar year and seconds from
                    # the beginning of that day.
                    # This converts Nagihara et al.'s time data into the format given
                    # in the original HFE dataset.

                    column_indicator = ".1" if probe == "p2" else ""
                    nagihara_data[mission][probe][sensor]["UTC_Time"] = np.vectorize(
                        days_since_year
                    )(
                        nagihara_spreadsheet[mission]["Day" + column_indicator],
                        int(mission[-4:]),
                    ) + np.vectorize(
                        seconds_interval
                    )(
                        nagihara_spreadsheet[mission]["Seconds" + column_indicator]
                    )
                    nagihara_data[mission][probe][sensor]["Time"] = np.vectorize(
                        a17_time
                    )(
                        np.vectorize(utc_to_tai)(
                            nagihara_data[mission][probe][sensor]["UTC_Time"]
                        )
                    )
                    # The original HFE dataset reduced temperature data into T and dT,
                    # where T is the average temperature across the bridge and dT is the
                    # difference between the upper and lower parts of the bridge.
                    # Nagihara et al. reduced the ALSEP data differently, explicitly
                    # giving temperature values for 'A' (upper) and 'B' (lower) parts of
                    # the bridge. This converts Nagihara et al.'s
                    # temperature data into the format given in the original HFE dataset.

                    TGA = nagihara_spreadsheet[mission][
                        "TG" + probe[1] + str(sensor) + "A"
                    ]
                    TGB = nagihara_spreadsheet[mission][
                        "TG" + probe[1] + str(sensor) + "B"
                    ]

                    # T

                    nagihara_data[mission][probe][sensor]["T"] = (TGA + TGB) / 2

                    # dT

                    nagihara_data[mission][probe][sensor]["dT"] = TGB - TGA

                    # We also retain Nagihara et al.'s explicitly-computed bridge values
                    # to avoid rounding errors later.

                    nagihara_data[mission][probe][sensor][
                        "TG" + probe[1] + str(sensor) + "A"
                    ] = TGA
                    nagihara_data[mission][probe][sensor][
                        "TG" + probe[1] + str(sensor) + "B"
                    ] = TGB

            # converts 'files' to pandas dataframes and drops empty rows.
            # we don't flag empty rows because they're just formatting artifacts.
            # previous versions wrote intermediate csv files, but we don't bother here;
            # may add later if efficiency becomes an issue or they prove desirable
            # for some other reason.

    for mission in nagihara_data:
        for probe in nagihara_data[mission]:
            for sensor in nagihara_data[mission][probe]:
                nagihara_data[mission][probe][sensor] = pd.DataFrame.from_dict(
                    nagihara_data[mission][probe][sensor]
                ).dropna()
    return nagihara_data


def ingest_nagihara_2019(nagihara_data=None, pathname="."):

    # grab files in directory according to Nagihara et al.'s naming convention
    # and read them as pandas dataframes

    # These files are parsable as csv with fields separated by a variable number of spaces.

    # as elsewhere, m is mission (here mission by year), p is probe number,
    # s is 'file,' i.e. bridge identifier.

    for mission in ["a15_1975", "a17_1975"]:
        nagihara_data[mission] = {}
        for probe in ["p1", "p2"]:
            nagihara_data[mission][probe] = {}
            for sensor in [1, 2]:
                filename = "{pathname}/{m}_hfe_1975_l2_arcsav_tg{p}{s}.tab".format(
                    pathname=pathname, m=mission[0:3], p=probe[1], s=sensor
                )
                if os.path.isfile(filename):
                    nagihara_data[mission][probe][sensor] = pd.read_csv(
                        filename, engine="python", sep=r"\ +"
                    )

                    # Time

                    # generate a separate column for mission epoch time.
                    # convert native DOY format (e.g. '1975-092T00:04:00.817')
                    # to datetime as an intermediate step.

                    # to_pydatetime is necessary because pandas otherwise gets mad
                    # about loss of (here nonexistent) nanosecond precision.

                    # then remove leap seconds for parity with NSSDC data and convert to mission epoch time.

                    nagihara_data[mission][probe][sensor]["UTC_Time"] = np.vectorize(
                        nagihara_doy_to_dt
                    )(nagihara_data[mission][probe][sensor]["time"])

                    if mission[0:3] == "a17":
                        nagihara_data[mission][probe][sensor]["Time"] = np.vectorize(
                            a17_time
                        )(
                            np.vectorize(tai_to_utc)(
                                nagihara_data[mission][probe][sensor][
                                    "UTC_Time"
                                ].dt.to_pydatetime()
                            )
                        )

                    elif mission[0:3] == "a15":
                        nagihara_data[mission][probe][sensor]["Time"] = np.vectorize(
                            a17_time
                        )(
                            np.vectorize(tai_to_utc)(
                                nagihara_data[mission][probe][sensor][
                                    "UTC_Time"
                                ].dt.to_pydatetime()
                            )
                        )
                    nagihara_data[mission][probe][sensor].drop(
                        columns=["time"], inplace=True
                    )

                    # dT

                    # Nagihara et al. include both the high- and low-sensitivity dT
                    # measurements for every data point. To construct our reduced sets,
                    # we choose one of these measurements per point as a canonical dT.
                    # For each point, we simply select the dT measurement that Nagihara
                    # et al. used to calculate their canonical bridge temperatures.

                    dT_init = nagihara_data[mission][probe][sensor]["Time"].size
                    nagihara_data[mission][probe][sensor]["dT"] = pd.Series(
                        np.zeros(dT_init)
                    )

                    TA = nagihara_data[mission][probe][sensor][
                        "TG{p}{s}A".format(p=probe[1], s=sensor)
                    ]
                    TB = nagihara_data[mission][probe][sensor][
                        "TG{p}{s}B".format(p=probe[1], s=sensor)
                    ]
                    DTH = nagihara_data[mission][probe][sensor][
                        "DTH{p}{s}".format(p=probe[1], s=sensor)
                    ]
                    DTL = nagihara_data[mission][probe][sensor][
                        "DTL{p}{s}".format(p=probe[1], s=sensor)
                    ]

                    # is TA - TB within a rounding error of a given measurement?
                    # then select that measurement as canonical dT, preferring DTH
                    # if both measurements qualify.

                    DTH_index = round(abs(TA - TB - DTH), 3) <= 0.01
                    DTL_index = round(abs(TA - TB - DTL), 3) <= 0.01

                    nagihara_data[mission][probe][sensor].loc[
                        DTL_index, "dT"
                    ] = DTL.loc[DTL_index]
                    nagihara_data[mission][probe][sensor].loc[
                        DTH_index, "dT"
                    ] = DTH.loc[DTH_index]

                    # are there points where neither measurement qualifies? something has
                    # gone wrong with the analysis.

                    if not np.all(np.bitwise_or(DTL_index, DTH_index)):
                        raise ValueError(
                            "Margin of error for PDS data dT selection appears to be off."
                        )

                    # then flip the sign for parity with the NSSDC convention

                    nagihara_data[mission][probe][sensor]["dT"] = (
                        nagihara_data[mission][probe][sensor]["dT"] * -1
                    )

                    # rename and reindex columns for parity with NSSDC data and later convenience

                    nagihara_data[mission][probe][sensor].rename(
                        columns={"TG{p}{s}avg".format(p=probe[1], s=sensor): "T"},
                        inplace=True,
                    )
                    nagihara_data[mission][probe][sensor] = nagihara_data[mission][
                        probe
                    ][sensor].reindex(
                        columns=[
                            "Time",
                            "T",
                            "dT",
                            "UTC_Time",
                            "TG{p}{s}A".format(p=probe[1], s=sensor),
                            "TG{p}{s}B".format(p=probe[1], s=sensor),
                        ]
                    )

    return nagihara_data


# functions for polishing dataset and splitting to thermometers.


def discard_flagged(data):
    data_clean = pd.DataFrame(columns=data.columns)
    # retain ambiguous dT corrections and time inversions
    forbidden_flags = sum(
        [
            0b1,
            0b100,
            0b1000,
            0b10000,
            0b100000,
            0b1000000,
            0b10000000,
            0b100000000,
            0b1000000000,
        ]
    )
    index = np.bitwise_and(data["flags"], forbidden_flags).values == 0
    for column in data.columns:
        data_clean[column] = data[column].loc[index].values
    return data_clean


def substitute_corrections(data):
    if "dT_corr" in data.columns:
        data["dT"] = data["dT_corr"]
        data.drop(columns=["dT_corr"], inplace=True)


def standardize_time(data, mission):

    # use directly-reformatted time from Nagihara; it's given to millisecond precision and
    # meaningful floating-point errors could plausibly be introduced.

    # otherwise go ahead and do the full conversion from epoch time.

    if "UTC_Time" in data.columns:
        data["Time"] = data["UTC_Time"]
        data.drop(columns=["UTC_Time"], inplace=True)
    elif mission[0:3] == "a17":
        data["Time"] = np.vectorize(a17_to_greg)(data["Time"])
        data["Time"] = np.vectorize(tai_to_utc)(data["Time"])
    elif mission[0:3] == "a15":
        data["Time"] = np.vectorize(a15_to_greg)(data["Time"])
        data["Time"] = np.vectorize(tai_to_utc)(data["Time"])


def polish_hfe(data):
    for mission in data.keys():
        for probe in data[mission].keys():
            for sensor in data[mission][probe].keys():
                data[mission][probe][sensor] = discard_flagged(
                    data[mission][probe][sensor]
                )
                substitute_corrections(data[mission][probe][sensor])
                standardize_time(data[mission][probe][sensor], mission)


def join_hfe_years(data):
    for mission in ["a15_1975", "a17_1975", "a17_1976", "a17_1977"]:
        for probe in data[mission]:
            for sensor in data[mission][probe]:
                data[mission[0:3]][probe][sensor] = (
                    data[mission[0:3]][probe][sensor]
                    .append(data[mission][probe][sensor])
                    .reset_index(drop=True)
                )
        del data[mission]


def split_to_thermometers(data):
    data_split = {}
    for mission in data.keys():
        data_split[mission] = {}
        for probe in data[mission].keys():
            data_split[mission][probe] = {}
            for sensor in data[mission][probe].keys():
                if mission[3:] == "":
                    # construct per-thermometer temperature fields for NSSDC data
                    probe_code = probe[1]
                    if sensor == 1:
                        sensor_code = ["G", "1"]
                    elif sensor == 2:
                        sensor_code = ["G", "2"]
                    elif sensor == 4:
                        sensor_code = ["R", "1"]
                    elif sensor == 5:
                        sensor_code = ["R", "2"]
                    if sensor != 3:
                        # calculate per-thermometer temperatures from NSSDC data
                        A = "T" + sensor_code[0] + probe_code + sensor_code[1] + "A"
                        B = "T" + sensor_code[0] + probe_code + sensor_code[1] + "B"
                        data_split[mission][probe][sensor] = pd.DataFrame(
                            columns=["Time", A, B]
                        )
                        data_split[mission][probe][sensor]["Time"] = data[mission][
                            probe
                        ][sensor]["Time"]
                        data_split[mission][probe][sensor][A] = (
                            data[mission][probe][sensor]["T"]
                            - 0.5 * data[mission][probe][sensor]["dT"]
                        )
                        data_split[mission][probe][sensor][B] = (
                            data[mission][probe][sensor]["T"]
                            + 0.5 * data[mission][probe][sensor]["dT"]
                        )
                        data_split[mission][probe][sensor]["flags"] = data[mission][
                            probe
                        ][sensor]["flags"]
                    elif sensor == 3:
                        TC = "TC{p}".format(p=probe_code)
                        data_split[mission][probe][sensor] = data[mission][probe][
                            sensor
                        ].drop("HTR", axis=1)
                        data_split[mission][probe][sensor] = data_split[mission][probe][
                            sensor
                        ].rename(
                            columns={
                                "TC1": TC + "1",
                                "TC2": TC + "2",
                                "TC3": TC + "3",
                                "TC4": TC + "4",
                                "TREF": "TREF",
                                "flags": "flags",
                            }
                        )
                    for column in data_split[mission][probe][sensor]:
                        # round time back to second precision rather than the
                        # microsecond precision introduced by astropy time scale conversion.
                        if column == "Time":
                            data_split[mission][probe][sensor][column] = data[mission][
                                probe
                            ][sensor][column].dt.round("1s")
                        # retains the (probably spurious) 5 digits after decimal given
                        # by the Lamont data, rather than the (definitely spurious)
                        # additional digits of precision introduced by Python's
                        # floating point representation.
                        elif column[0] == "T":
                            data_split[mission][probe][sensor][column] = round(
                                data_split[mission][probe][sensor][column], 5
                            )
                else:
                    # simply use temperatures calculated by Nagihara et al. to avoid introducing errors
                    data_split[mission][probe][sensor] = data[mission][probe][
                        sensor
                    ].drop(columns=["T", "dT"])
                    for column in data_split[mission][probe][sensor]:
                        # retain millikelvin precision of Nagihara et al. sets rather than
                        # spurious precision introduced by numpy floating point representation
                        if column[0:2] == "TG":
                            data_split[mission][probe][sensor][column] = round(
                                data_split[mission][probe][sensor][column], 3
                            )
    return data_split


# functions & depth dictionary for writing combined 'depth' set

# cm below surface for each sensor. 'S' marks sensors lying somewhere on the surface,
# 'N' marks the fourth thermocouple at the same position as the top gradient thermometer,
# 'X' marks off-scale bridges. the assembly functions for the third / depth set
# exclude these ranges from consideration; they are included for possible future work.

DEPTHDICT = {
    "a17": {
        "p1": {
            1: {"TG11A": 130, "TG11B": 177},
            2: {"TG12A": 185, "TG12B": 233},
            3: {"TC11": "S", "TC12": 14, "TC13": 66, "TC14": "N", "TREF": "S"},
            4: {"TR11A": 140, "TR11B": 167},
            5: {"TR12A": 195, "TR12B": 223},
        },
        "p2": {
            1: {"TG21A": 131, "TG21B": 178},
            2: {"TG22A": 186, "TG22B": 234},
            3: {"TC21": "S", "TC22": 15, "TC23": 67, "TC24": "N", "TREF": "S"},
            4: {"TR21A": 140, "TR21B": 169},
            5: {"TR22A": 196, "TR22B": 224},
        },
    },
    "a15": {
        "p1": {
            1: {"TG11A": 35, "TG11B": 84},
            2: {"TG12A": 91, "TG12B": 139},
            3: {"TC11": "S", "TC12": "S", "TC13": 0, "TC14": "N", "TREF": "S"},
            4: {"TR11A": 45, "TR11B": 73},
            5: {"TR12A": 101, "TR12B": 129},
        },
        "p2": {
            1: {"TG21A": [-6, "X"], "TG21B": [32, "X"]},
            2: {"TG22A": 49, "TG22B": 97},
            3: {"TC21": "S", "TC22": "S", "TC23": "S", "TC24": 0, "TREF": "S"},
            4: {"TR21A": [3, "X"], "TR21B": [41, "X"]},
            5: {"TR22A": 59, "TR22B": 87},
        },
    },
}

# concatenates all subsurface sensors per probe (barring TC4) and assigns depth values.


def combine_with_depth(data):
    data_clean = {}
    for mission in ["a15", "a17"]:
        data_clean[mission] = {}
        for probe in ["p1", "p2"]:
            data_clean[mission][probe] = pd.DataFrame(
                columns=["Time", "T", "sensor", "depth", "flags"]
            )
    for mission in data.keys():
        for probe in data[mission].keys():
            for sensor in data[mission][probe].keys():
                for column in data[mission][probe][sensor].columns:
                    # is it a temperature value from a sensor we want to include in this set?
                    if column[0] == "T" and column[1] != "i":
                        if (
                            isinstance(DEPTHDICT[mission][probe][sensor][column], int)
                            and DEPTHDICT[mission][probe][sensor][column] > 0
                        ):
                            depth_slice = data[mission][probe][sensor][
                                ["Time", column, "flags"]
                            ]
                            depth_slice = depth_slice.reindex(
                                ["Time", "T", column, "sensor", "depth", "flags"],
                                axis=1,
                            )
                            depth_slice["depth"] = DEPTHDICT[mission][probe][sensor][
                                column
                            ]
                            depth_slice["sensor"] = column
                            depth_slice["T"] = depth_slice[column]
                            depth_slice = depth_slice.drop(column, axis=1)
                            data_clean[mission][probe] = data_clean[mission][
                                probe
                            ].append(depth_slice)
            data_clean[mission][probe] = (
                data_clean[mission][probe].sort_values(by="Time").reset_index(drop=True)
            )
    return data_clean
