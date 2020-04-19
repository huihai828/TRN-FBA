class ParserError(Exception):
    def __init__(self, line, msg):
        self.nline = str(line)
        self.msg = msg
        #self.line = line

# class myError(Exception):
#     pass
    # def __init__(self, nline, msg):
    #     #self.line = line