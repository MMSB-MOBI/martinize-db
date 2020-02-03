import Express from 'express';

export enum ErrorType {
  /** SERVER ERRORS */
  Server = 1, 

  /** CLIENT ERRORS: Not Found */
  NotFound = 101, 

  /** CLIENT ERRORS: Unauthorized access */
  Forbidden = 201,
  TokenInvalid,

  /** CLIENT ERRORS: Bad Request */
  Format = 301,
  MissingParameters,
}

const ErrorsToText = {
  [ErrorType.Server]: [500, ""],
  [ErrorType.NotFound]: [404, ""],
  [ErrorType.Forbidden]: [403, ""],
  [ErrorType.TokenInvalid]: [403, ""],
  [ErrorType.Format]: [400, ""],
  [ErrorType.MissingParameters]: [400, ""],
};

export default new class Errors {
  send(code: ErrorType, res: Express.Response) {
    const [http_code, message] = ErrorsToText[code] as [number, string];

    res.status(http_code).json({
      code,
      message
    });
  }

  throw(code: ErrorType) : never {
    const [http_code, message] = ErrorsToText[code] as [number, string];

    throw new ApiError(String(http_code), {
      code,
      message
    });
  }
};

export class ApiError extends Error {
  constructor(public message: string, public data: any) { 
    super(message); 
    this.name = "ApiError";
  }
}
