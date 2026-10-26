from ionotec.download import ibge
from ionotec.download import cors
from ionotec.download import geonet
from ionotec.download import sugar
from ionotec.download import sonel
from ionotec.download import hongkong
from ionotec.download import igs
from ionotec.download import garner
from ionotec.download import earthscope
from ionotec.download import australia
from ionotec.download import ramsac


from pathlib import Path

import pandas as pd
import logging

logging.basicConfig(
    filename="sources_errors.log",
    level=logging.ERROR,
    format="%(asctime)s %(message)s",
)


def getSource(station):
    df_station_source = pd.read_csv('stations_source.csv').set_index('station')
    return df_station_source.loc[station,'source']



class downrinex:

    def __init__(self, start, end, output_folder='.'):

        self.date_start = start
        self.date_end = end
        self.dict_station_source = {}
        self.list_stations = []

        self.folder_deposit = output_folder

        Path(output_folder).mkdir(parents=True, exist_ok=True)


    def list_available_stations(self,source='all'):
        self.dict_station_source = {}

        #### IGS: RINEX 4
        if source in ['all','igs']:
            try:
                available_stations = igs.stations_available_igs(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='igs'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "igs: stations_available_igs failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING igs unavailable, skipped (see sources_errors.log)")

        #### EARTHSCOPE: RINEX 4
        if source in ['all','earthscope']:
            try:
                available_stations = earthscope.stations_available_earthscope(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='earthscope'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "earthscope: stations_available_earthscope failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING earthscope unavailable, skipped (see sources_errors.log)")

        #### AUSTRALIA: RINEX 3
        if source in ['all','australia']:
            try:
                available_stations = australia.stations_available_australia(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='australia'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "australia: stations_available_garner failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING australia unavailable, skipped (see sources_errors.log)")

        #### GARNER: RINEX 3
        if source in ['all','garner']:
            try:
                available_stations = garner.stations_available_garner(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='garner'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "garner: stations_available_garner failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING garner unavailable, skipped (see sources_errors.log)")


        #### IBGE (Brasil): RINEX 3
        if source in ['all', 'ibge']:
            try:
                available_stations = ibge.stations_available_ibge(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station] = 'ibge'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "ibge: stations_available_ibge failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING ibge unavailable, skipped (see sources_errors.log)")

        #### GEONET (New Zealand): RINEX 3
        if source in ['all','geonet']:
            try:
                available_stations = geonet.stations_available_geonet(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='geonet'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "geonet: stations_available_geonet failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING geonet unavailable, skipped (see sources_errors.log)")

        #### HONGKONG (Sumatra): RINEX 3
        if source in ['all','hongkong']:
            try:
                available_stations = hongkong.stations_available_hongkong(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='hongkong'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "hongkon: stations_available_hongkong failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING hongkong unavailable, skipped (see sources_errors.log)")


        #### CORS: RINEX 2
        if source in ['all', 'cors']:
            try:
                available_stations = cors.stations_available_cors(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station] = 'cors'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "cors: stations_available_cors failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING cors unavailable, skipped (see sources_errors.log)")


        #### SUGAR: RINEX 2
        if source in ['all','sugar']:
            try:
                available_stations = sugar.stations_available_sugar(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='sugar'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "sugar: stations_available_sugar failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING sugar unavailable, skipped (see sources_errors.log)")

        #### SONEL: RINEX ??
        if source in ['all','sonel']:
            try:
                available_stations = sonel.stations_available_sonel(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='sonel'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "sonel: stations_available_sonel failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING sugar unavailable, skipped (see sources_errors.log)")

        #### ARGENTINA: RINEX 2
        if source in ['all','ramsac']:
            try:
                available_stations = ramsac.stations_available_ramsac(self.date_start, self.date_end)
                for station in available_stations:
                    if station in self.list_stations: continue
                    self.dict_station_source[station]='ramsac'
                    self.list_stations.append(station)
            except Exception:
                logging.exception(
                    "ramsac: stations_available_ramsac failed (%s - %s)",
                    self.date_start, self.date_end,
                )
                print("WARNING ramsac unavailable, skipped (see sources_errors.log)")


        df = pd.DataFrame(self.dict_station_source.items(), columns=['station', 'source'])
        df.to_csv('stations_source.csv')

        return self.list_stations

    def download_rinex_station(self, station):
        source = self.dict_station_source[station]
        if source  == 'ibge':
            return ibge.download_ibge(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'cors':
            return cors.download_cors(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'geonet':
            return geonet.download_geonet(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'sugar':
            return sugar.download_sugar(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'sonel':
            return sonel.download_sonel(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'hongkong':
            return hongkong.download_hongkong(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'igs':
            return igs.download_igs(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'garner':
            return garner.download_garner(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'earthscope':
            return earthscope.download_earthscope(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'australia':
            return australia.download_australia(station, self.date_start, self.date_end, output_dir=self.folder_deposit)
        if source  == 'ramsac':
            return ramsac.download_ramsac(station, self.date_start, self.date_end, output_dir=self.folder_deposit)



