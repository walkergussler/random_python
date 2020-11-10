#!/usr/bin/python
from SimpleHTTPServer import SimpleHTTPRequestHandler
import BaseHTTPServer
import pdb
import json
import dnds

class CORSRequestHandler (SimpleHTTPRequestHandler):
    # def end_headers (self):
    #     self.send_header('Access-Control-Allow-Origin', '*')
    #     self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
    #     self.send_header('Access-Control-Allow-Headers', '*')
    #     SimpleHTTPRequestHandler.end_headers(self)

    def _set_headers(self):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', '*')
        self.end_headers()

    def do_OPTIONS(self):           
        self.send_response(200, "ok")       
        self.send_header('Access-Control-Allow-Origin', '*')                
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With")        

    def do_GET(self):           
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Content-type',    'text/html')                                    
        self.end_headers()              
        self.wfile.write("<html><body>Hello world!</body></html>")
        self.connection.shutdown(1) 

    def do_POST(self):
        content_length = int(self.headers['Content-Length']) # <--- Gets the size of data
        post_data = self.rfile.read(content_length) # <--- Gets the data itself

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("STARTING: sliding window, gaussian-smoothed dN/dS analysis...")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        qry_and_ref_seq_list = post_data.replace("s1=","").replace("s2=","")
        qry_and_ref_seq_list = qry_and_ref_seq_list.split("&")

        qry_seq_raw = qry_and_ref_seq_list[0]
        ref_seq_raw = qry_and_ref_seq_list[1]

        dnds_data = dnds.dnds_pipeline(qry_seq_raw, ref_seq_raw)

        dnds_data_vec = dnds_data[0]
        qry_seq_indices = dnds_data[1]
        
        #return_data = '{"dnds_data_vec": dnds_data_vec, "qry_seq_indices": qry_seq_indices}'

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("ENDING: sliding window, gaussian-smoothed dN/dS analysis...")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        self._set_headers()
        self.wfile.write("<html><body><h1>"+str(dnds_data_vec)+"</h1><h2>"+str(qry_seq_indices)+"</h2></body></html>")
        #self.wfile.write(str(dnds_data))


        

if __name__ == '__main__':
    BaseHTTPServer.test(CORSRequestHandler, BaseHTTPServer.HTTPServer)