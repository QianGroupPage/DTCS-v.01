from datetime import datetime, timezone
import json
from typing import List, Dict

import pymongo
from bson.objectid import ObjectId
import monty
import pandas as pd
import numpy as np

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
    
    def store(self, xps: XPSExperiment, ts: CRNTimeSeries, fake=False) -> Dict:
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

        if not fake:
            self._mongo[crn_collection].insert_one(doc)

        return doc

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
        input_env = xps_spectra.to_numpy()
        for i, doc in enumerate(docs):
            data = CRNData(doc)
            db_env = np.array([])
            raw_db_env = data.xps_data.envelope.to_numpy()
            input_i = 0

            scale_factor = max(input_env) / max(raw_db_env)
            raw_db_env *= scale_factor

            # Resample the stored data to match the number of datapoints in the external
            for i, pt in enumerate(raw_db_env):
                if pt > input_env[input_i]:
                    if i == 0:
                        db_env = np.append(db_env, pt)
                    else:
                        db_env = np.append(db_env, raw_db_env[i-1])
                    input_i += 1
                if input_i == len(input_env):
                    break

            if len(db_env) < len(input_env):
                db_env = np.append(db_env, np.array([db_env[-1] for _ in range(len(input_env)-len(db_env))]))

            rmse = ((input_env - db_env)**2).mean() **.5
            if rmse < min_rmse or min_rmse < 0:
                min_rmse = rmse
                min_rmse_data = data

        print("Min RMSE:", min_rmse)
        return min_rmse_data
