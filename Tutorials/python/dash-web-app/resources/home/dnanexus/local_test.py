# -*- coding: utf-8 -*-
import my_app

app = my_app.create_app()

if __name__ == '__main__':
    app.run_server(debug=True)

