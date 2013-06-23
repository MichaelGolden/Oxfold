package uk.ac.ox.osscb.util;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.filter.Filter;
import ch.qos.logback.core.spi.FilterReply;

public class StdOutFilter extends Filter<ILoggingEvent>{

	@Override
	public FilterReply decide(ILoggingEvent logEvt) {		
		if(logEvt.getLevel().isGreaterOrEqual(Level.WARN)){
			return FilterReply.DENY;
		}
		return FilterReply.ACCEPT;			
	}
}
