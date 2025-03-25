
########################################################################
class Container:
    def __init__(self,title=""):
        self.title = title
        self.container = np.array([], dtype=object)
        self.titles = {}
        self.size = 0
    
    def __add__(self,other):
        if not isinstance(other,Collection):
            raise TypeError("unsupported types of operands!")
        return self.copy().extend(other.get())
        
    def __sub__(self,other):
        if not isinstance(other,Collection):
            raise TypeError("unsupported types of operands!")
        oCopy = self.copy()
        for key in other.get_titles():
            if key in list(self.titles.keys()):
                oCopy.__delitem__(key)
        return oCopy
        
    def __len__(self):
        return self.size
        
    def __contains__(self,key):
        if isinstance(key,str):
            if not self.has(key):
                return False 
        elif isinstance(key,int):
            if key > self.size:
                return False
        else:
            raise ValueError("Container must be refered either by text title or integer index!")
        return True
    
    def __iter__(self):
        if not len(self.container):
            return iter([])
        #return iter(list(map(lambda obj: obj.Obj, self.container)))
        return iter([obj.Obj for obj in self.container])
    
    def __getitem__(self,index):
        if isinstance(index, slice):
            return self.get(index.start,index.stop)
        return self.get(index)
    
    def __setitem__(self,key,Obj):
        index,title = self._parse_key(key)
        if title=="~attribute":
            setattr(self, key, Obj)
            return
        link = None
        if index:
            link = self.container[index-1].Obj.link
        Obj = ContainerElement(Obj,index,link)
        self.titles[title] = Obj
        self.container[index] = Obj
    
    def __delitem__(self,key):
        index,title = self._parse_key(key)
        self.container = np.delete(self.container,index)
        del self.titles[title]
        self.size -= 1
        
    def _parse_key(self,key):
        if not key:
            key = 0
        if isinstance(key,str):
            if not self.has(key):
                if hasattr(self,key):
                    return key,"~attribute"
                raise ValueError(f"There is no object with title {key}!") 
            index = self.titles[key].index
            title = key
        elif isinstance(key,int):
            index = key
            title = self.container[key].title
        else:
            raise ValueError("Container must be refered either by text title or integer index!")
        return index,title
    
    def _get_title(self,obj):
        if hasattr(obj, 'title'):
            return obj.title
        return ""

    def _check_object_title(self,Obj):
        if not hasattr(Obj,'title') or not Obj.title:
            Obj.title = str(self.size)
        if self.has(Obj.title):
            raise KeyError(f"Title {Obj.title} is already occupied!")
        return True
            
    def clear(self):
        self.container = np.array([], dtype=object)
        self.titles = {}
        self.size = 0

    def has(self,title):
        return title in list(self.titles.keys())
    
    def index(self,title):
        if not self.has(title):
            return -1
        return self.titles[title].Obj.index
    
    def append(self,Obj):
        self._check_object_title(Obj)
        link = None
        if self.size:
            link = self.container[self.size-1].link
        Obj = ContainerElement(Obj,self.size,link)
        self.titles[Obj.title] = Obj
        self.container = np.append(self.container, Obj)
        self.size += 1
        
    def insert(self,key,Obj):
        index,title = self._parse_key(key)
        if abs(index) >= self.size:
            self.container.append(Obj)
            return
        self._check_object_title(Obj)
        Obj = ContainerElement(Obj,index,None)
        if index == 0:
            self.container[0].link = Obj.link_function
        else:
            Obj.link = self.container[index-1].link_function
            self.container[index].link = Obj.link_function
        self.container = np.insert(self.container, index, Obj)
        self.titles[Obj.title] = Obj
        self.size += 1            
        
    def extend(self,ls):
        for Obj in ls:
            if hasattr(Obj,'title') and self.has(Obj.title):
                continue
            self.append(Obj)
        
    def get_titles(self):
        return list(self.titles.keys())
    
    def get(self,start_key=None,stop_key=None):
        if start_key==None and stop_key==None:
            #return list(map(lambda obj: obj.Obj.copy(), self.container))
            return [obj.Obj.copy() for obj in self.container]
        index_start,title_start = self._parse_key(start_key)
        if title_start=="~attribute":
            return eval(compile(f"self.{index_start}", "<string>", "eval"))
        if stop_key:
            index_stop,title_stop = self._parse_key(stop_key)
            indices = [index_start,index_stop]
            indices.sort()
            #return list(map(lambda i: self.container[i].Obj, range(indices[0],indices[1]+1,1)))
            return [self.container[i].Obj for i in range(indices[0],indices[1]+1,1)]
        else:
            return self.container[index_start].Obj
            
    def dict(self):
        return self.titles
        
    def push(self,obj_ls):
        if self.size and obj_ls and type(self[0]) != type(obj_ls[0]):
            raise TypeError("Object in the collection and the list are of different types!")
        self.clear()
        for Obj in obj_ls:
            self.append(Obj)
            
    def sort(self,key="",reverse=False):
        obj_list = self.get()
        self.clear()
        # if key==None, sorting is escaped
        if key != None:
            if key=="":
                obj_list.sort(key=lambda Obj: Obj.title)
            else:
                try:
                    eval(compile(f"obj_list.sort(key={key})", "<string>", "eval"))
                except:
                    raise TypeError(f"command {key} cannot be used for sorting these objects!")
        if reverse:
            obj_list.reverse()
        for Obj in obj_list:
            self.append(Obj)
            
    def sorted(self,key="",reverse=False):
        return self.copy().sort(key,reverse)
        
    def reverse(self):
        return self.sorted(key=None,reverse=True)
            
    def copy(self):
        oNewContainer = Collection(self.title)
        for i in range(self.size):
            try:
                oNewContainer.append(self.container[i].Obj.copy())
            except:
                oNewContainer.append(copy.deepcopy(self.container[i].Obj))
        return oNewContainer

########################################################################
class Collection:
    def __init__(self, title=""):
        self.title = title
        self.container = []
        self.para = {}
        
    def _get_key(self,key):
        if type(key)==type(0):
            if len(self) <= key:
                return
            return key
        elif type(key)==type(""):
            return self.index(key)
    
    def _get_title(self,obj):
        if hasattr(obj, 'title'):
            return obj.title
        return ""

    def __len__(self):
        return len(self.container)
        
    def __contains__(self,key):
        if self._get_key(key) != None:
            return True
        return False
    
    def __iter__(self):
        if not self.container:
            return iter([])
        records = []
        for record in self.container:
            records.append(record)
        return iter(records)
    
    def __getitem__(self,key):
        key = self._get_key(key)
        if key != None:
            return self.container[key]
    
    def __setitem__(self,key,value):
        key = self._get_key(key)
        if key != None:
            self.container[key] = value
    
    def __delitem__(self,key):
        key = self._get_key(key)
        if key != None:
            del self.container[key]
            
    def __repr__(self):
        #return "\n".join(list(map(lambda i: "%d\t%s" % (i+1,str(self.container[i])), range(len(self.container)))))
        return "\n".join([f"{i + 1}\t{self.container[i]}" for i in range(len(self.container))])
            
    def __str__(self):
        #return ";".join(list(map(lambda Obj: str(Obj), self.container)))
        return ";".join([str(Obj) for Obj in self.container])
            
    def has(self,title):
        return title in self.get_titles()
    
    def index(self,title):
        if type(title) == type(0):
            index = title
            if abs(index) >= len(self):
                return
            return index
        titles = self.get_titles()
        if title in titles:
            return titles.index(title)
        return
    
    def append(self,obj):
        self.container.append(obj)
        
    def extend(self,ls):
        self.container.extend(ls)
        
    def get_titles(self):
        #return list(map(lambda obj: obj.title, self.container))
        return [obj.title for obj in self.container]
    
    def get(self,titles=[]):
        if titles:
            #return list(filter(lambda Obj: Obj.title in titles, self.container))
            return [Obj for Obj in self.container if Obj.title in titles]
        else:
            return self.container
            
    def copy(self):
        oCollection = Collection(self.title)
        oCollection.para.update(self.para)
        for record in self.container:
            oCollection.append(record.copy())
        return oCollection
            
########################################################################
class ContainerElement:
    def __init__(self,Obj,index,link=None):
        if not hasattr(Obj,'title'):
            Obj.title = str(index)
        self.title = Obj.title
        self.index = index
        self.link = link
        self.Obj = Obj
        
    def link_function(self,command,args=[]):
        if command == "increment":
            self.index += 1
        elif command == "decrement":
            self.index -= 1
        else:
            pass
            
