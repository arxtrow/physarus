use std::{backtrace::Backtrace, error::Error, fmt, sync::Arc};

#[macro_export]
macro_rules! pbrt_err {
    // Basic error
    ($msg:expr) => {
        $crate::core::PbrtError::new($msg)
    };

    // With formatting
    ($fmt:expr, $($arg:tt)*) => {
        $crate::core::PbrtError::new(format!($fmt, $($arg)*))
    };
}

#[derive(Debug)]
pub struct PbrtError {
    pub message: String,
    pub source: Option<Arc<dyn Error + Send + Sync>>,
    pub contexts: Vec<String>,
    pub backtrace: Backtrace,
}

impl PbrtError {
    pub fn new<M: Into<String>>(message: M) -> Self {
        Self {
            message: message.into(),
            source: None,
            contexts: vec![],
            backtrace: Backtrace::capture(),
        }
    }

    pub fn with_source<E>(mut self, err: E) -> Self
    where
        E: Error + Send + Sync + 'static,
    {
        self.source = Some(Arc::new(err));
        self
    }

    pub fn context<C: Into<String>>(mut self, ctx: C) -> Self {
        self.contexts.push(ctx.into());

        self
    }
}

impl fmt::Display for PbrtError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "\n\n*******************ERRRRORRR********************")?;
        writeln!(f, "{}", self.message)?;

        if self.contexts.len() > 0 {
            writeln!(f, "\nCaused by:")?;
            for (i, ctx) in self.contexts[0..].iter().enumerate() {
                writeln!(f, "    {}: {}", i, ctx)?;
            }
        }

        writeln!(f, "************************************************\n\n")?;

        Ok(())
    }
}

impl Error for PbrtError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.source.as_ref().map(|arc| arc.as_ref() as _)
    }
}
