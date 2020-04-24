#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import dxpy
import dash
import dash_core_components as dcc
import dash_html_components as html
import my_app

@dxpy.entry_point('main')
def main():

    app = my_app.create_app()
    dxhandler = dxpy.get_handler(dxpy.JOB_ID)
    dxhandler.set_properties({"httpsAppState": "running"})
    app.run_server(host='0.0.0.0', port=443)

    return 1

dxpy.run()
