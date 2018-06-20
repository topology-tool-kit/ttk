#ifndef _PDBARYCENTERIMPL_H
#define _PDBARYCENTERIMPL_H


template <class dataType> int ttk::PersistenceDiagramsBarycenter::execute() const{

	Timer t;
	{

	for(int i=0; i<numberOfInputs_; i++){
		std::cout << "#" <<std::endl;
		std::vector<diagramTuple>* CTDiagram = static_cast<std::vector<diagramTuple>*>(inputData_[i]);
		for (unsigned int j=0; j<CTDiagram->size(); j++)
		{
			diagramTuple t = CTDiagram->at(j);
			dataType persistence = std::get<4>(t);
			std::cout << persistence <<std::endl;
		}	
	}
	
		
	std::stringstream msg;
	msg << "[PersistenceDiagramsBarycenter] processed in "
		<< t.getElapsedTime() << " s. (" << threadNumber_
		<< " thread(s))."
		<< std::endl;
	dMsg(std::cout, msg.str(), timeMsg);
	}

	return 0;
}

#endif
