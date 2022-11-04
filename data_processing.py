import numpy as np

def report(title, ans, units):
    try:
        print (title, ": ("+ans.report()+")", units)
    except:
        print (title,":", ans[0], units, "±", ans[1], units)

def weighted_average(stuff):
    # tups is list of tuples: val, err as outputted by max model
    try:
        vals = [result.val for result in stuff]
        weights = 1/np.array([result.tot for result in stuff]) ** 2
        #print (vals)
        #print (weights)
    except:
        vals = [tup[0] for tup in stuff]
        weights = 1/np.array([tup[1] for tup in stuff]) ** 2
    return np.dot(vals, weights) / np.sum(weights), 1/np.sqrt(np.sum(weights))

def sum(*args):
    sum = 0
    errs2 = 0
    try:
        for tup in args:
            sum += tup[0]
            errs2 += tup[1]**2
    except:
        for tup in args[0]:
            sum += tup[0]
            errs2 += tup[1]**2
    return sum, np.sqrt(errs2)

def minus(tup):
    return -tup[0], tup[1]

def resultify(tup):
    return Result(tup[0], stat = tup[1])

class Result:

    def __init__(self, val, stat = 0, sys = 0, tot = 0):
        self.val = val
        self.stat = stat
        self.sys = sys
        self.tot = np.sqrt(stat**2 + sys**2) if tot == 0 else tot 

    def report(self):
        string = str(self.val)
        if self.stat != 0:
            string += " ± " + str(self.stat) + " (stat)"
        if self.sys != 0:
            string += " ± " + str(self.sys) + " (sys)"
        if self.stat == 0 and self.sys == 0 and self.tot != 0:
            string += " ± " + str(self.tot) + " (total)"
        return string

    def __str__(self):
        return self.report()

    def sin(self):
        val = np.sin(self.val)
        if self.stat == 0 and self.sys == 0:
            tot = abs(np.cos(self.val) * self.tot)
            return Result(val, tot = tot)
        stat = abs(np.cos(self.val) * self.stat)
        sys = abs(np.cos(self.val) * self.sys)
        return Result(val, stat = stat, sys = sys)

    def cos(self):
        return (self*-1 + np.pi / 2).sin()
        
    
    def __add__(self, other):
        try:
            val = self.val + other.val
            stat = np.sqrt(self.stat**2 + other.stat**2)
            sys = np.sqrt(self.sys**2 + other.stat**2)  
        except:
            val = self.val + other
            stat = self.stat
            sys = self.sys
        tot = np.sqrt(stat**2 + sys**2)
        return Result(val, stat = stat, sys = sys, tot = tot)

    def __mul__(self, other):
        try:
            val = self.val * other.val
            stat = val * np.sqrt((self.stat/self.val)**2 + (other.stat/other.val)**2)
            sys = val * np.sqrt((self.sys/self.val)**2 + (other.sys/other.val)**2)
        except:
            val = self.val * other
            stat = self.stat * other
            sys = self.sys * other
        tot = np.sqrt(stat**2 + sys**2)
        return Result(val, stat = stat, sys = sys, tot = tot)   
    
    def __sub__(self, other):
        return self + other*(-1)

    def __truediv__(self, other):
        try:
            val = self.val / other.val
            stat = val * np.sqrt((self.stat/self.val)**2 + (other.stat/other.val)**2)
            sys = val * np.sqrt((self.sys/self.val)**2 + (other.sys/other.val)**2)
        except:
            val = self.val / other
            stat = self.stat / other
            sys = self.sys / other
        tot = np.sqrt(stat**2 + sys**2)
        return Result(val, stat = stat, sys = sys, tot = tot)  

    def __pow__(self, other):
        try:
            val = self.val ** other.val
            stat = np.sqrt((other.val * self.val**(other - 1) * self.stat)**2 + (self.val**other.val * np.log(self.val) * other.stat)**2)
            sys = np.sqrt((other.val * self.val**(other - 1) * self.sys)**2 + (self.val**other.val * np.log(self.val) * other.sys)**2)
        except:
            val = self.val ** other
            stat = other * self.val**(other-1) * self.stat
            sys = other * self.val**(other-1) * self.stat
        tot = np.sqrt(stat**2 + sys**2)
        return Result(val, stat = stat, sys = sys, tot = tot)

    def __gt__(self, other):
        try:
            return self.val > other.val
        except:
            return self.val > other

    def __le__(self, other):
        return not self > other

    def __lt__(self, other):
        try:
            return self.val < other.val
        except:
            return self.val < other
    
    def __ge__(self, other):
        return not self < other

    def __eq__(self, other):
        try:
            return self.val == other.val
        except:
            return self.val == other

    def equals(self, other):
        try:
            return ((self <= other) and (self.val + self.tot > other.val - other.tot)) or ((self > other) and (self.val - self.tot < other.val + other.tot))
        except:
            return ((self <= other) and (self.val + self.tot > other)) or ((self > other) and (self.val - self.tot < other))
