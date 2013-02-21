#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# gatk_realign_and_recalibrate 0.0.1 test suite
# Generated by dx-app-wizard.
#
# See http://wiki.dnanexus.com/Developer-Tutorials/Build-Your-First-DNAnexus-App
# for instructions on how to modify this file.

import os, sys, unittest, json, subprocess

import dxpy, dxpy.app_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")

def makeInputs():
    # Please fill in this method to generate default inputs for your app.
    return {}

class Testgatk_realign_and_recalibrate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_input = makeInputs()
        bundled_resources = dxpy.app_builder.upload_resources(src_dir)
        cls.app_id = dxpy.app_builder.upload_app(src_dir, bundled_resources, overwrite=True)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_base_input(self):
        job = dxpy.DXApp(self.app_id).run(self.base_input)
        print "Waiting for job to complete"
        job.wait_on_done()
        print json.dumps(job.describe()["output"])
