from datetime import datetime, timezone
import json
from typing import List

import pymongo
from bson.objectid import ObjectId
import monty
import pandas as pd

from lblcrn.experiments.xps import XPSExperiment
from lblcrn.experiments.time_series import CRNTimeSeries

crn_collection = "crnsims"

class CRNData:
    def __init__(self, doc) -> None:
        self.xps_data = pd.read_json(doc["xps_data"],convert_axes=False)
        self.ts_data = pd.read_json(doc["ts_data"], convert_axes=False)
        self.rsys = doc["rsys"]
        self.rsys_id = doc["rsys_id"]

class CRNStorage:
    def __init__(self, uri: str="mongodb://localhost:27017", db: str="lblcrn") -> None:
        self._mongo = pymongo.MongoClient(uri)[db]
    
    def store(self, xps: XPSExperiment, ts: CRNTimeSeries) -> None:
        raw_ts = ts.df.to_json(double_precision=15)
        raw_xps = xps.df.simulated.to_json(double_precision=15)
        rsys_fingerprint = ts.rsys.fingerprint()
        rsys_id = ts.rsys.id()

        doc_id = ObjectId()
        timestamp = datetime.now(timezone.utc)
        doc = {
            "_id": doc_id,
            "rsys_id": rsys_id,
            "rsys_fingerprint": rsys_fingerprint,
            "rsys": str(ts.rsys),
            "xps_data": raw_xps,
            "ts_data": raw_ts,
            "created_at": timestamp,
        }

        self._mongo[crn_collection].insert_one(doc)

    def load(self, doc_id: str) -> CRNData:
        doc = self._mongo[crn_collection].find_one({"_id": ObjectId(doc_id)})
        return CRNData(doc)

    def load_from_rsys_fingerprint(self, rsys) -> List[CRNData]:
        docs = self._mongo[crn_collection].find({"rsys_fingerprint": rsys.fingerprint()})
        
        data = []
        for doc in docs:
            data.append(CRNData(doc))

        return data

    def load_from_rsys_id(self, rsys) -> List[CRNData]:
        docs = self._mongo[crn_collection].find({"rsys_id": rsys.id()})
        
        data = []
        for doc in docs:
            data.append(CRNData(doc))

        return data


    def find_closest_xps(self, rsys, xps_spectra: pd.Series) -> CRNData:
        """Given the rsys and xps spectra, find and return the closest match in the database.

        The closest match is defined by the minimum RMSE between the spectra.
        """
        docs = self._mongo[crn_collection].find({"rsys_fingerprint": rsys.fingerprint()})
        
        min_rmse = -1
        min_rmse_data = None
        for i, doc in enumerate(docs):
            data = CRNData(doc)
            rmse = ((xps_spectra.to_numpy() - data.xps_data.envelope.to_numpy())**2).mean() **.5
            if rmse < min_rmse or min_rmse < 0:
                min_rmse = rmse
                min_rmse_data = data

        print("Min RMSE:", min_rmse)
        return min_rmse_data
