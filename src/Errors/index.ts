import Express from 'express';

export enum ErrorType {
  /** SERVER ERRORS */
  Server = 1, 

  /** CLIENT ERRORS: Not Found */
  NotFound = 101, 
  UserNotFound,
  ElementNotFound,

  /** CLIENT ERRORS: Unauthorized access */
  Forbidden = 201,
  TokenInvalid,
  InvalidPassword,
  Unallowed,

  /** CLIENT ERRORS: Bad Request */
  Format = 301,
  MissingParameters,
  UsernameExists,
  EmailExists,
  InvalidMethod,
  TooManyFiles,
  FileTooLarge,
  InvalidMoleculeFiles,
}

const ErrorsToText = {
  [ErrorType.Server]: [500, "Server error"],
  [ErrorType.NotFound]: [404, "Page not found"],
  [ErrorType.UserNotFound]: [404, "User not found"],
  [ErrorType.ElementNotFound]: [404, "Element not found"],
  [ErrorType.Forbidden]: [403, "Access restricted"],
  [ErrorType.Unallowed]: [403, "You don't have the right to do that"],
  [ErrorType.TokenInvalid]: [403, "Invalid or expired token"],
  [ErrorType.InvalidPassword]: [403, "Invalid password"],
  [ErrorType.Format]: [400, "Parameter format is invalid"],
  [ErrorType.MissingParameters]: [400, "Missing parameters"],
  [ErrorType.TooManyFiles]: [400, "Too many files specified"],
  [ErrorType.FileTooLarge]: [400, "Sended file is too large"],
  [ErrorType.InvalidMoleculeFiles]: [400, "Sended molecule files (ITP, PDB/GRO) are incorrect"],
  [ErrorType.UsernameExists]: [409, "Username already exists"],
  [ErrorType.EmailExists]: [409, "Email already exists"],
  [ErrorType.InvalidMethod]: [405, "Method not allowed"],
};

export default new class Errors {
  send(code: ErrorType, res: Express.Response, additionnal?: object) {
    if (res.headersSent) {
      return this.throw(code);
    }

    const [http_code, message] = ErrorsToText[code] as [number, string];

    res.status(http_code).json({
      code,
      message,
      ...(additionnal || {})
    });
  }

  /**
   * Throw the given error.
   * If you're in a promise, make sure the error is correctly catched and sent !
   */
  throw(code: ErrorType, additionnal?: object) : never {
    const [http_code, message] = ErrorsToText[code] as [number, string];

    throw new ApiError(String(http_code), {
      code,
      message,
      ...(additionnal || {})
    });
  }

  /**
   * Make the given error.
   * If you're in a promise, make sure the error is correctly catched and sent !
   */
  make(code: ErrorType, additionnal?: object) : ApiError {
    const [http_code, message] = ErrorsToText[code] as [number, string];

    return new ApiError(String(http_code), {
      code,
      message,
      ...(additionnal || {})
    });
  }
};

export class ApiError extends Error {
  constructor(public message: string, public data: any) { 
    super(message); 
    this.name = "ApiError";
  }
}
