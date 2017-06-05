#!/usr/bin/python

import sys, copy, tarfile

""" - Splits Percolator output into decoy and target files.
    - Extracts unique PSM/peptides/proteins out of a Percolator output file.
    - Merges Percolator output files
    Usage: python percolator_output_modifier.py command psm/peptides/proteins [score] infile outfile [outfile2]
    """

try:
    from lxml import etree
except Exception:
    sys.stderr.write('Failed to import lxml module.')


def readPercout(fname):
    doc = None
    try:
        doc = etree.parse(fname)
    except Exception:
        sys.stderr.write('Could not parse XML provided in %s or error reading file. \n' % (fname))
    return doc


def splitTargetDecoy(doc, args, ns):
    """ Splits XML into target/decoy/notspecified elements.
        then calls refill function to create a new XML ElementTree with results.
        
        Usage: splitTargetDecoy('test.xml', ['psms', 'peptides', 'proteins'])
    """
    for arg in args:
        assert arg in ['psms', 'peptides', 'proteins'], Exception('Filtering argument must be one or more of: psms, peptides, proteins.')

    output_elts = {
    'psms' : doc.xpath('//xmlns:psms', namespaces=ns),
    'peptides' : doc.xpath('//xmlns:peptides', namespaces=ns),
    'proteins' : doc.xpath('//xmlns:proteins', namespaces=ns)
    }

    def do_split(name, tree, ns):
        """ Does the actual target/decoy splitting. Returns lists of classified elements."""
        # check if not passed an empty tree
        if len(tree) == 0:
            sys.stdout.write('Cannot output %s since it is not in the provided input XML. Continuing.\n' % name)
            return [etree.Element('none')]*3
        
        # check if decoy attribute is in tree at all
        if tree[0].xpath('//@xmlns:decoy', namespaces=ns) == []:
            sys.stderr.write('Percolator output has no specified target/decoy spectra, unable to split.\n')
            return [etree.Element('none')]*3
        
        # filter
        ttree = []
        dtree = []
        discarded = []
#        ttree = etree.Element(name, ns)
#        dtree = etree.Element(name, ns)
#        discarded = etree.Element(name, ns)
        searchElt = name[:-1]
        
        for element in tree[0].xpath('xmlns:%s' % searchElt, namespaces=ns):
            if element.xpath('@xmlns:decoy', namespaces=ns) == ['true']:
                dtree.append(element)
            elif element.xpath('@xmlns:decoy', namespaces=ns) == ['false']:
                ttree.append(element)
            else:
                discarded.append(element)

        return ttree, dtree, discarded
    
    # call splitfunction, then check output
    target,decoy,discarded = {},{},{}
    for arg in args:
        target[arg], decoy[arg], discarded[arg] = do_split(arg, output_elts[arg], ns)
    for arg in args:
        if [target[arg], decoy[arg], discarded[arg]] != 3*['none']:
            break
        sys.stderr.write('No data matching %s has been found in Percolator output \
                            file.\n' % [x for x in args])

    # clean parsed elementtree and make new ones
    target = refillTree(doc, output_elts.keys(), target, ns)
    decoy = refillTree(doc, output_elts.keys(), decoy, ns)
    discarded = refillTree(doc, output_elts.keys(), discarded, ns)
   
#    sys.stdout.write('%s\n' % etree.tostring(target))

    return target, decoy #, discarded?


def refillTree(doc, oldelements, newelements, ns):
    """Takes an ElementTree, takes out specified oldelements (by tag). Replaces
    with list of newelements. Returns a new ElementTree.
    """
    # remove specified old elements
    newdoc = copy.deepcopy(doc)
    root = newdoc.getroot()
    elements_toremove = []
    for el in oldelements:
        removes = root.xpath('xmlns:%s' % el, namespaces=ns) # element 'psms' or 'peptides'
        for remove in removes:
            children = remove.getchildren()
            for el in children:
                remove.remove(el)
    
    # put in new ones
    for node in root.getchildren():
        try:
            for child in newelements[node.tag[node.tag.index('}')+1:]]:
                node.append(child)
        except: # PSMs, peptides should be in this order and not in newelements-dict's arbitrary order.
            pass
    return etree.ElementTree(root)



def filterUniques(tar, to_filter, score, ns):
    """ Filters unique psms/peptides/proteins from (multiple) Percolator output XML files.
        Takes a tarred set of XML files, a filtering query (e.g. psms), a score to
        filter on and a namespace.
        Outputs an ElementTree.
    """
    for tf in to_filter:
        assert tf in ['psms', 'peptides', 'proteins'], Exception('filterUnique function needs a specified to_filter list of psms, peptides, proteins.')
    assert score in ['q','pep','p'], Exception('filterUnique function needs a specified score to filter on of q, pep or p.')
    try:
        with tarfile.open(tar, 'r') as f:
            members = f.getmembers()
            f.extractall()
    except:
        sys.stderr.write('Could not extract Percolator files from dataset: %s \n' % tar)
        return 1
    
    docs = []
    for fn in members:
        docs.append(etree.parse(fn.name))

    # lookup dicts
    scores = {'q':'q_value', 'pep':'pep', 'p':'p_value'} 
    filt_el_dict = {'psms':'xmlns:peptide_seq', 'peptides':'@xmlns:peptide_id' }

    # result dict
    filtered = {'psms':{}, 'peptides':{}, 'proteins':{} }
    for doc in docs:
        for filt_el in to_filter:
            feattree = doc.xpath('//xmlns:%s' % filt_el, namespaces=ns)
            if feattree == []:
                sys.stdout.write('%s not found in (one of the) Percolator output documents. Continuing...\n' % filt_el)
                continue
            for feat in feattree[0]:
                # It's actually faster to loop through the feat's children, 
                # but this is 2-line code and still readable.
                featscore = float(feat.xpath('xmlns:%s' % scores[score], namespaces=ns)[0].text)
                seq = feat.xpath('%s' % filt_el_dict[filt_el], namespaces=ns)
                
                try: # psm seqs are parsed here
                    seq = seq[0].attrib['seq']
                except Exception: ## caught when parsing peptide seqs (different format)
                    seq = str(seq[0])
                
                if seq not in filtered[filt_el]:
                    filtered[filt_el][seq] = feat
                elif featscore < filtered[filt_el][seq]:
                    #FIXME now it only works for LOWER than scores (eg q-vals, pep, but not for scores that are better when higher)
                    filtered[filt_el][seq] = feat
        
    # make trees from filtered dicts
    for filt_el in filtered:
        outlist = []
        for feat in filtered[filt_el].values():
            outlist.append(feat)
        filtered[filt_el] = outlist
#        node = etree.Element(filt_el)
#        node.extend(filtered[filt_el].values())
#        filtered[filt_el] = node
    
    outdoc = refillTree(docs[0], ['psms', 'peptides', 'proteins'], filtered, ns)

    return outdoc


def mergeXML(datasets, ns):
    tomerge = ['psms', 'peptides', 'proteins']
    tomerge_el = ['psms', 'peptides', 'proteins']
    outdoc = readPercout(datasets[0])
    root = outdoc.getroot()
    for el in tomerge:
        tomerge_el.extend(root.xpath('//xmlns:%s' % el, namespaces=ns))
    for fn in datasets[1:]:
        doc = readPercout(fn)
        for tm, tme in zip(tomerge, tomerge_el):
            print tm,tme,"both"
            for el in doc.xpath('//xmlns:%s' % tm, namespaces=ns):
                print el,"el"
                tme.append(el)
            print tme,"tme"
    print tomerge_el,"to merge el"
    return outdoc

def untar(tar):
    try:
        with tarfile.open(tar, 'r') as f:
            members = f.getmembers()
            f.extractall()
            return [x.name for x in members]
    except:
        sys.stderr.write('Could not extract Percolator files from dataset: %s \n' % tar)
        return 1

def writeXML(*args):
    """Takes a filename and _ElementTree arguments in tuples ((fn, tree), (fn2,tree2),...) 
    and writes an XML file.
    """
    
    try:
    # hello! shouldn't these be exceptions instead?
        assert False not in [isinstance(x, tuple) for x in args] # arguments are tuples
        assert False not in [len(x)>2 for x in args] # tuples are 2-length
        assert False not in [isinstance(x[0], str) for x in args] # first arg from tuple is string for fn.
        assert False not in [isinstance(x[1], etree._ElementTree) for x in args] # second arg from tuple is ElementTree.
    except:
        Exception('function writeXML takes arguments in form of tuples (filename, XMLtree)')
        
    for fn, doc in args:
        try:
            outfile = open(fn, 'w')
        except Exception:
            sys.stderr('Unable to write XML output to file %s \n' % fn)
        doc.write(outfile)


def parseOptions(args):
    ns = {'xmlns':'http://per-colator.com/percolator_out/13', 'xmlns:p' : 'http://per-colator.com/percolator_out/13', \
    'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance', \
    'xsi:schemaLocation':'http://per-colator.com/percolator_out/13,https://github.com/percolator/percolator/raw/pout-1-3/src/xml/percolator_out.xsd', \
    'p:majorVersion':'2', 'p:minorVersion':'04', 'p:percolator_version':'Percolator version 2.04'}
    
    if args[0] == 'splittd':
        fname = args[2]
        doc = readPercout(fname)
        target, decoy = splitTargetDecoy(doc, args[1].split(','), ns)
        writeXML((args[3], target), (args[4], decoy))
    
    elif args[0] == 'filter_uni':
        tosplit = args[2].split(',')
        uniques = filterUniques(args[1], tosplit, args[3], ns)
        writeXML((args[4], uniques))
    
    elif args[0] == 'merge': # is there any use for this?
        with open(args[1]) as fp:
            datasets = fp.read()
        print datasets,"these are the datasets"
        datasets = [opt.strip() for opt in datasets.strip().split('\n')]
        sys.stdout.write('\ndatasets: %s' % datasets)
        merged = mergeXML(datasets, ns)
        writeXML((args[2], merged))
    
    elif args[0] == 'untar_merge':
        files = untar(args[1])
        sys.stdout.write('%s'%files)
        merged = mergeXML(files, ns)
        writeXML((args[2], merged))
        
    else:
        sys.stderr.write('Argument %s not recognized.\n' % args[0])
        return 1


def main():
    parseOptions(sys.argv[1:])


if __name__ == '__main__':
    main()
