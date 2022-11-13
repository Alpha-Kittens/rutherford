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
        return np.dot(vals, weights) / np.sum(weights), 1/np.sqrt(np.sum(weights))
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

def resultify_fit(model, params):
    pass

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
            sys = np.sqrt(self.sys**2 + other.sys**2)  
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
            #print (self)
            #print (other)
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
            stat = np.sqrt((other.val * self.val**(other.val - 1) * self.stat)**2 + (self.val**other.val * np.log(self.val) * other.stat)**2)
            sys = np.sqrt((other.val * self.val**(other.val - 1) * self.sys)**2 + (self.val**other.val * np.log(self.val) * other.sys)**2)
        except:
            val = self.val ** other
            stat = abs(other) * self.val**(other-1) * self.stat
            sys = abs(other) * self.val**(other-1) * self.stat
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

def attempts(options):
    if len(options) == 0:
        return [[]]
    else:
        ret = []
        prior = attempts(options[1:])
        if type(options[0]) != type("hello"):
            for option in options[0]:
                for attempt in prior:
                    ret.append([option] + attempt)
        else:
            for attempt in prior:
                ret.append([options[0]] + attempt)
        return ret

class AsymmetricResult(Result):
    def __init__(self, val, stathigh = 0, statlow = 0, syshigh = 0, syslow = 0):
        super().__init__(val, stat = max(stathigh, statlow), sys = max(syshigh, syslow))
        self.stathigh = stathigh
        self.statlow = statlow
        self.syshigh = syshigh
        self.syslow = syslow

    def report(self):
        string = str(self.val)
        if self.stathigh != 0 or self.statlow != 0:
            if self.stathigh == self.statlow and self.stathigh != 0:
                string += " ± " + str(self.stat) + " (stat)"
            else:
                string += " +" + str(self.stathigh) + "/-" + str(self.statlow) + " (stat)"
            
        if self.syshigh != 0 or self.syslow != 0:
            if self.syshigh == self.syslow and self.syshigh != 0:
                string += " ± " + str(self.sys) + " (sys)"
            else:
                string += " + " + str(self.syshigh) + "/- " + str(self.syslow) + " (sys)"

        return string

    @staticmethod
    def asymmetric_evaluate(function, *args, progress = False):
        """
        Given a function and arguments, some of which are Result objects, finds error bounds on the evaluation of the function.
        Assumes that critical points of function are on boundary of uncertainty region; reasonable for well-behaved functions and smallish relative errors?
        Arguments:
            * `function`: Function to be evaluated
            * `*args`: Arguments to be passed to `function`, in order. At least one should be a `Result`. Otherwise, what are you even doing lol
        Returns:
            * `result`: AsymmetricResult object. Value is evaluation of `function` at `*args`; error bounds are calculated by testing `funciton` on error bounds of `Results` in `*args`.
        NOTE: AsymmetricResult is a subclass of Result, and can be treated as one. AsymmetricResults store the upper and lower stat/sys errors separately; the corresponding
            `Result` object uses the max(lower stat/sys error, upper stat/sys error) as its error for its stat/sys uncertainty. 
        """
        highstat = 0
        lowstat = 0
        highsys = 0
        lowsys = 0
        stats = []
        syss = []
        eval_args = []
        for arg in args:
            try:
                stats.append((arg.val + arg.stat, arg.val - arg.stat))
                syss.append((arg.val + arg.sys, arg.val - arg.sys))
                eval_args.append(arg.val)
            except:
                stats.append((arg))
                syss.append((arg))
                eval_args.append(arg)
        val = function(*eval_args)
        stat_attempts = attempts(stats)
        sys_attempts = attempts(syss)
        total_attempts = len(stat_attempts) + len(sys_attempts)
        i = 1
        for stat_attempt in stat_attempts:
            if progress:
                print (str(i)+"/"+str(total_attempts))
                i += 1
            try:
                eval = function(*stat_attempt)
                if eval > val + highstat:
                    highstat = abs(eval - val)
                if eval < val - lowstat:
                    lowstat = abs(eval - val)
            except Exception as e:
                print (e, ", skipping")

        for sys_attempt in sys_attempts:
            if progress:
                print (str(i)+"/"+str(total_attempts))
                i += 1
            try:
                eval = function(*sys_attempt)
                if eval > val + highsys:
                    #print ("highsys")
                    highsys = abs(eval - val)
                if eval < val - lowsys:
                    #print ("lowsys")
                    lowsys = abs(eval - val)
            except Exception as e:
                print (e, ", skipping")
        return AsymmetricResult(val, stathigh = highstat, statlow = lowstat, syshigh = highsys, syslow = lowsys)


