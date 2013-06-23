package uk.ac.ox.osscb.util;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.filter.Filter;
import ch.qos.logback.core.spi.FilterReply;
// 
	
/**
 * {@see <a href="http://logback.qos.ch/manual/filters.html">http://logback.qos.ch/manual/filters.html</a>}
 * @author Vladimir
 *
 */
public class StdErrFilter  extends Filter<ILoggingEvent>{

	@Override
	public FilterReply decide(ILoggingEvent logEvt) {
		if(logEvt.getLevel().isGreaterOrEqual(Level.WARN)){
			return FilterReply.ACCEPT;			
		}
		return FilterReply.DENY;
	}
}
