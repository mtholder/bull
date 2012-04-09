#!/usr/bin/env python

import sys, os, random, math
from cStringIO import StringIO
_RNG = random.Random()

_logger_initialized = False

_user_cfg_dir = os.path.expanduser("~/.dendropy")

def get_dendropy_cfg_file(file_name):
    """Looks for file_name in _user_cfg_dir, if it is not found there then
    the package directory is checked.
    
    Returns `None` if the file is not found.
    """
    global _user_cfg_dir
    if _user_cfg_dir:
        fp = os.path.join(_user_cfg_dir, file_name)
        if os.path.exists(fp):
            return True, fp
    from pkg_resources import resource_filename
    fp = resource_filename(__name__, file_name)
    if os.path.exists(fp):
        return True, fp
    return False, fp

def get_logger(s):
    """Wrapper around logging.getLogger that make sure that the dendropy
    logging configuration file is read (or a default is applied)
    """
    global _logger_initialized
    import logging
    msg = ""
    if not _logger_initialized:
        import logging.config
        filename = "logging_simtree.conf"
        found, full_path = get_dendropy_cfg_file(filename)
        specified_path_failed = False
        if found and os.path.exists(full_path):
            try:
                logging.config.fileConfig(full_path)
                _logger_initialized = True
                logger = logging.getLogger(s)
                return logger
            except:
                specified_path_failed = True
        else:
            msg = "logging config file not found at %s" % full_path
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        default_fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(default_fmt_str)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        if specified_path_failed:
            warning = 'Could not parse %s, using default logging' % full_path
            sys.stderr.write(warning)
        _logger_initialized = True
    l = logging.getLogger(s)
    if msg:
        l.debug(msg)
    return l
_LOG = get_logger("simtree")


class Node(object):
    def __init__(self, parent=None, rate_of_evol=1.0):
        self.name = str(id(self))
        self.par = parent
        self.children = []
        self.duration = 0.0
        self.ini_rate = rate_of_evol
        self.final_rate = None
        self.rate = self.ini_rate

    def finalize_rate(self, roeroe, max_rate, min_rate):
        global _RNG
        assert self.final_rate is None
        lnstddev = roeroe*self.duration
        if lnstddev > 0.0:
            # setting mean to ini_rate - stddev/2
            mu = math.log(self.ini_rate) - (lnstddev/2.0)
            _LOG.debug("Ln rate from LogN(%g, %g), In rate %g" % (mu, lnstddev, self.ini_rate))
            r = _RNG.lognormvariate(mu, lnstddev)
            self.final_rate =  min(max_rate, max(r, min_rate))
            _LOG.debug("Ln Rate of %g (orig %g) from LogN(%g, %g). Rate %g -> %g" % (math.log(self.final_rate), math.log(r), mu, lnstddev, self.ini_rate, self.final_rate))
        else:
            self.final_rate = self.ini_rate
        self.rate = (self.final_rate + self.ini_rate)/2.0

    def speciate(self, roeroe, max_rate=1e300, min_rate=1e-300):
        self.finalize_rate(roeroe, max_rate, min_rate)
        self.children.append(Node(parent=self, rate_of_evol=self.final_rate))
        self.children.append(Node(parent=self, rate_of_evol=self.final_rate))
        return self.children

    def prune_self(self):
        if self.par:
            p = self.par
            p.children.remove(self)
            if not p.children:
                p.prune_self()

    def add_time(self, t):
        self.duration += t

    def suppress_deg_two(self):
        self.rate = (self.final_rate + self.ini_rate)/2.0
        n = 0
        for c in self.children:
            n += c.suppress_deg_two()
        if len(self.children) == 1:
            c = self.children[0]
            sd = self.duration + c.duration
            wt = self.duration*self.rate + c.duration*c.rate
            self.rate = wt/sd
            self.duration = sd
            self.final_rate = c.final_rate
            self.children = c.children
            for gc in self.children:
                gc.par = self
            n += 1
            if c.name:
                self.name = c.name
        return n

    def write_newick(self, o, brlens=True):
        if self.children:
            f = True
            o.write("(")
            for c in self.children:
                if f:
                    f = False
                else:
                    o.write(",")
                c.write_newick(o)
            o.write(")")
        elif self.name:
            o.write(self.name)
        if brlens and self.par:
            o.write(":%f" % (self.duration*self.rate))

    def __str__(self):
        return self.newick(brlens=False)

    def newick(self, brlens=True):
        o = StringIO()
        self.write_newick(o, brlens=brlens)
        return o.getvalue()
    
def remove_random(l):
    n = len(l)
    ind = _RNG.randrange(n)
    el = l[ind]
    l.pop(ind)
    return el
    
  
def sim_tree(leaves_to_gen, options):
    birth = options.birth
    if birth <= 0.0:
        sys.exit("Speciation (birth) rate must be positive")
    death = max(0.0, options.death)
    irate = options.ini
    if irate <= 0.0:
        sys.exit("The initial rate of evolution must be positive")
    roeroe = options.rate
    if roeroe < 0.0:
        sys.exit("--rate value cannot be negative")
    eventrate = birth + death
    dprob = death/eventrate
    min_rate = options.min_rate
    if min_rate <= 0.0:
        sys.exit("The minimum rate of evolution must be positive")
    max_rate = options.max_rate
    if max_rate < min_rate:
        sys.exit("The maximum rate of evolution must be greater than the minimum rate.")

    root = Node(rate_of_evol=irate)
    evolving = list(root.speciate(roeroe, max_rate, min_rate))
    _LOG.debug(root.newick(brlens=False))
    ntips = 2
    while ntips < leaves_to_gen:
        total_rate = ntips*eventrate
        waiting_time = _RNG.expovariate(total_rate)
        for i in evolving:
            i.add_time(waiting_time)
        to_mod = remove_random(evolving)
        if _RNG.random() < dprob:
            to_mod.prune_self()
            _LOG.debug("Deleting %s after %g" % (to_mod.name, waiting_time))
            _LOG.debug(root.newick(brlens=False))
            ntips -= 1
            if ntips == 0:
                return None, []
        else:
            _LOG.debug("Speciating %s after %g" % (to_mod.name, waiting_time))
            evolving.extend(to_mod.speciate(roeroe, max_rate, min_rate))
            _LOG.debug(root.newick(brlens=False))
            ntips += 1
    # add a duration to all terminal branches.  Duration is U[0,1] * waiting time to the next event
    total_rate = ntips*eventrate
    next_ev = _RNG.expovariate(total_rate)
    curr_time = _RNG.random() * next_ev
    for i in evolving:
        i.add_time(curr_time)
        i.finalize_rate(roeroe, max_rate, min_rate)
    return root, evolving
            
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=0, 
        type="int",
        help="The random number generator seed")
    parser.add_option("-n", "--num", dest="n", default=1, 
        type="int",
        help="The number of trees to generate")
    parser.add_option("-l", "--leaves", dest="ntax", default=4, 
        type="int",
        help="The number of leaves to generate")
    parser.add_option("-c", "--sample", dest="nchosen", default=-1, 
        type="int",
        help="The number of leaves to sample (-1 to sample all leaves generated)")
    parser.add_option("-b", "--birth", dest="birth", default=1.0, 
        type="float",
        help="The birth (speciation rate). default 1.")
    parser.add_option("-d", "--death", dest="death", default=0.0, 
        type="float",
        help="The death (extinction rate). default 0.")
    parser.add_option("-r", "--rate", dest="rate", default=0.0, 
        type="float",
        help="This number will be multiplied by the branch length to determine the standard deviation of the log of the rate of evolution. The rate for the branch will then be the mean of the values at the endpoints.")
    parser.add_option("-i", "--initial-rate", dest="ini", default=1.0, 
        type="float",
        help="The initial rate of evolution for the root of the tree.")
    parser.add_option("-f", "--min-rate", dest="min_rate", default=1e-300, 
        type="float",
        help="The minimum rate of evolution")
    parser.add_option("-x", "--max-rate", dest="max_rate", default=1e300, 
        type="float",
        help="The maximum rate of evolution")
    (options, args) = parser.parse_args()
    if args:
        sys.exit("No arguments accepted.  Use the -h flag to see the list of command line options (flags) that are accepted")
    if options.seed < 1:
        import time
        options.seed = int(time.time()*1000)
    _LOG.debug("%s -s%d -n%d -l%d -c%d -b%g -d%g -r%g -i%g -f%g -x%f\n" % (
                     sys.argv[0], 
                     options.seed,
                     options.n,
                     options.ntax,
                     options.nchosen,
                     options.birth,
                     options.death,
                     options.rate,
                     options.ini,
                     options.min_rate,
                     options.max_rate
                     ))
    _RNG.seed(options.seed)
    n = options.n
    leaves_to_gen = options.ntax
    leaves_to_include = options.nchosen
    if leaves_to_include < 0:
        leaves_to_include = leaves_to_gen
    if leaves_to_include < 2:
        if leaves_to_include == 0:        
            print "();"
        else:
            print "(1);"
        sys.exit(0)
    if leaves_to_include > leaves_to_gen:
        sys.exit("--sample value cannot be larger than --leaves value.")
    for i in xrange(n):
        r = None
        while r is None:
            r, tips = sim_tree(leaves_to_gen, options)
        rleaves = _RNG.sample(tips, leaves_to_include)
        if leaves_to_include < leaves_to_gen:
            for l in tips:
                if not (l in rleaves):
                    _LOG.debug("Deleting %s" % l.name)
                    l.prune_self()
                    _LOG.debug(r.newick(brlens=False))
        for n, l in enumerate(rleaves):
            l.name = "%d" % (n+1)
        r.suppress_deg_two()
        print str(r)
